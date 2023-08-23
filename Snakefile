## Pipeline dependencies - samtools, bcftools, guppy, and whatshap
import pandas as pd

## For ftp files
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

## Pipeline Config
configfile: "config.yaml"

## Reference Genomes
REFIDS=['hs37d5','GRCh38']


## Loading flowcell metadata
flowcells = pd.read_csv("flowcell_metadata.csv").set_index("flowcell_id")
GIABIDS=set(list(flowcells["giab_id"]))

## Running pipeline
rule all:
    input: 
        expand("qc/{flowcell}_{refid}.bam.stats.tsv.gz", 
               refid = REFIDS,
               flowcell = list(flowcells.index)),
        expand("phased_{refid}.bam", refid = REFIDS), 
         expand("phased_{refid}.bam.bai", refid = REFIDS)

def get_config(wildcards):
    return flowcells.loc[wildcards.flowcell,"config"]

def get_fq(wildcards):
	fq_path = flowcells.loc[wildcards.flowcell, 'fq']
	return(fq_path)

## QC fastq files
rule fastqc:
    input: get_fq
    output:
        html="qc/fastqc/{flowcell}.html",
        zip="qc/fastqc/{flowcell}.zip"
    wrapper:
        "0.27.1/bio/fastqc"

## Aligning to reference -----------------------------------------------------------
def get_readgroup(wildcards):
    ## Format bam read group 
    run      = flowcells.loc[wildcards.flowcell, 'sample_id']
    model    = flowcells.loc[wildcards.flowcell, 'platform']
    lib      = flowcells.loc[wildcards.flowcell, 'sample_id']
    date     = flowcells.loc[wildcards.flowcell, 'date']
    gpy_cfg  = flowcells.loc[wildcards.flowcell,'config']
    sample   = flowcells.loc[wildcards.flowcell, 'giab_id']
    
    read_group =   (
        f"@RG\\tID:{run}\\t"
        f"PU:{wildcards.flowcell}\\t"
        f"PL:nanopore\\t"
        f"PM:{model}\\t"
        f"LB:{lib}\\t"
        f"DT:{date}\\t"
        f"PG:guppy-{config['guppy_version']}-{gpy_cfg}\\t"
        f"SM:{sample}"
    )
    
    return(read_group)

rule get_hs37d5:
    input: FTP.remote(config["hs37d5"]["ref"])
    output: "resources/hs37d5.fna"
    shell: "gunzip -c {input} > {output}"

rule get_GRCh38:
    input: FTP.remote(config["GRCh38"]["ref"])
    output: "resources/GRCh38.fna"
    shell: "gunzip -c {input} > {output}"    

rule index_ref:
    input: "resources/{refid}.fna"
    output: "resources/{refid}.fna.fai"
    wrapper: "0.38.0/bio/samtools/faidx"

rule map_reads:
    input:
        ref="resources/{refid}.fna", 
        refidx = "resources/{refid}.fna.fai",
        fastq=get_fq
    output: "bams/{flowcell}_{refid}.bam"
    params: read_group=get_readgroup, mem=48, threads=22
    threads: 22
    conda: "envs/map_reads.yaml"
    shell: """
        minimap2 -t {params.threads} -aL -z 600,200 -x map-ont \
                -R \'{params.read_group}\' {input.ref} {input.fastq} \
            | samtools sort -m {params.mem}G -@{params.threads} \
                -O bam --reference {input.ref} > {output}
        samtools index {output}
    """
    
## BAM QC 
rule samtools_stats:
    input:
        "bams/{flowcell}_{refid}.bam"
    output:
        "samtools_stats_{refid}/{flowcell}_{refid}.txt"
    log:
        "logs/samtools_stats/{flowcell}_{refid}.log"
    wrapper:
        "0.38.0/bio/samtools/stats"

rule bam_stats:
    input:
        "bams/{flowcell}_{refid}.bam"
    output:
        "qc/{flowcell}_{refid}.bam.stats.tsv.gz"
    conda: "envs/bam_stats.yaml"
    shell:
        "python scripts/get_bam_stat.py {input} {output}"
    
rule multiqc:
    input:
        expand("samtools_stats_{{refid}}/{flowcell}_{{refid}}.txt",
               flowcell = list(flowcells.index))
    output:
        "qc/multiqc_{refid}.html"
    log:
        "logs/multiqc_{refid}.log"
    wrapper:
        "0.38.0/bio/multiqc"
## Phasing reads --------------------------------------------------------------------
def get_phased_vcf_url(wildcards):
	vcf_url= config[wildcards.ref_id]["vars"][wildcards.giab_id]
	return(vcf_url)

rule get_phased_vcf:
	output: 
		vcf="resources/{ref_id}_{giab_id}.vcf.gz",
		vcf_idx="resources/{ref_id}_{giab_id}.vcf.gz.tbi",
	params: vcf_url=get_phased_vcf_url
	shell: """
		curl -o {output.vcf} {params.vcf_url}
		tabix -p vcf {output.vcf}
	"""

rule phase_bams:
    input: 
        bam = "bams/{flowcell}_{refid}.bam", 
        bamidx = "bams/{flowcell}_{refid}.bam",
        ref = "resources/{refid}.fna",
        refidx = "resources/{refid}.fna.fai",
        var = "resources/{refid}.vcf.gz",
        varidx = "resources/{refid}.vcf.gz.tbi",
    output: "bams/{flowcell}_{refid}_phased.bam"
    log: "logs/phase_bams_{flowcell}_{refid}.log"
    conda: "envs/phase.yaml"
    shell: """
        whatshap haplotag \
            -o {output} \
	         -r {input.ref} \
            --ignore-read-groups \
	        --skip-missing-contigs \
            {input.var} \
            {input.bam} &> {log}
    """

rule index_phased:
    input: rules.phase_bams.output
    output: "bams/{flowcell}_phased_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"


## Combine bams --------------------------------------------------------------------
def get_bams(wildcards):
    bams = expand("bams/{flowcell}_{{refid}}_phased.bam", flowcell = list(flowcells.index))
    return(bams)

rule combine_bams:
    input: get_bams
    output: "phased_{refid}.bam"
    log: "logs/combine_bams_{refid}.log"
    params: "-f -O bam"
    threads: 8
    wrapper: "0.38.0/bio/samtools/merge"

rule index_combined:
    input: "phased_{refid}.bam"
    output: "phased_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"
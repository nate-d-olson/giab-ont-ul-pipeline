## Pipeline dependencies - samtools, bcftools, guppy, and whatshap
import pandas as pd

## For ftp files
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

## Pipeline Config
configfile: "config.yaml"

## Reference Genomes
REFIDS=['hs37d5','GRCh38']
#REFIDS=['GRCh38']


## Loading flowcell metadata
flowcells = pd.read_csv("flowcell_metadata.csv").set_index("flowcell_id")
GIABIDS=set(list(flowcells["giab_id"]))

## Running pipeline
rule all:
    input: 
        expand("qc/{giabid}_{flowcell}_GRCh38.bam.stats.tsv.gz", 
               zip,
               giabid = list(flowcells["giab_id"]),
               flowcell = list(flowcells.index)),
        expand("qc/{giabid}_{flowcell}_hs37d5.bam.stats.tsv.gz", 
               zip,
               giabid = list(flowcells["giab_id"]),
               flowcell = list(flowcells.index)),
         expand("{giabid}_phased_{refid}.bam", giabid=GIABIDS, refid = REFIDS), 
         expand("{giabid}_phased_{refid}.bam.bai", giabid=GIABIDS, refid = REFIDS)
        ## Multi QC not recognizing output
        # expand("qc/multiqc_{refid}.html", refid = REFIDS)

        
def get_config(wildcards):
    return flowcells.loc[wildcards.flowcell,"config"]

## Basecalling --------------------------------------------------------------------
## Skipping basecalling for UCSC data
# rule basecalling:
#     input: "../../raw/{flowcell}"
#     output: directory("fastq/{flowcell}")
#     params: get_config
#     shell: 'sw/ont-guppy/bin/guppy_basecaller -x cuda:0 -r \
#                 --num_callers 8 \
#                 --gpu_runners_per_device 4 \
#                 --chunks_per_runner 1664 \
#                 --input_path {input} \
#                 --save_path {output} \
#                 --records_per_fastq 0 \
#                 --config {params}'

# rule combine_and_compress_fastq:
#     input: rules.basecalling.output
#     output: "fastq/{flowcell}.fastq.gz"
#     shell: 'cat {input}/*fastq | gzip > {output}'

# ## QC fastq files
# rule fastqc:
#     input: "fastq/{flowcell}.fastq.gz"
#     output:
#         html="qc/fastqc/{flowcell}.html",
#         zip="qc/fastqc/{flowcell}.zip"
#     wrapper:
#         "0.27.1/bio/fastqc"

## Aligning to reference -----------------------------------------------------------
def get_readgroup(wildcards):
    ## Format bam read group 
    run      = flowcells.loc[wildcards.flowcell, 'sample_id']
    model    = flowcells.loc[wildcards.flowcell, 'platform']
    lib      = flowcells.loc[wildcards.flowcell, 'sample_id']
    date     = flowcells.loc[wildcards.flowcell, 'date']
    gpy_cfg  = flowcells.loc[wildcards.flowcell,'config']
    flowcell = flowcells.loc[wildcards.flowcell, 'type']
    kit      = flowcells.loc[wildcards.flowcell, 'kit']
    sample   = flowcells.loc[wildcards.flowcell, 'giab_id']
    
    read_group =   (
        f"@RG\\tID:{run}\\t"
        f"PU:{wildcards.flowcell}\\t"
        f"PL:nanopore\\t"
        f"PM:{model}\\t"
        f"LB:{lib}\\t"
        f"DT:{date}\\t"
        f"PG:guppy-{config['guppy_version']}-{gpy_cfg}\\t"
        f"DS:Flowcell={flowcell},kit={kit}\\t"
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

def get_fq(wildcards):
	fq_path = flowcells.loc[wildcards.flowcell, 'fq']
	return(fq_path)

rule map_reads:
    input:
        ref="resources/{refid}.fna", 
        refidx = "resources/{refid}.fna.fai",
        fastq=get_fq
    output: "bams/{giab_id}_{flowcell}_{refid}.bam"
    threads: 7
    log: "logs/map_reads_{giab_id}_{flowcell}_{refid}.log"
    params: read_group=get_readgroup, mem=96
    conda: "envs/map_reads.yaml"
    shell: """
        minimap2 -t {threads} -aL -z 600,200 -x map-ont \
                -R \'{params.read_group}\' {input.ref} {input.fastq} \
            | samtools sort -m {params.mem}G \
                -O bam --reference {input.ref} > {output} 2>{log}
        samtools index {output} &>>{log}
    """
    
## BAM QC 
rule samtools_stats:
    input:
        "bams/{giab_id}_{flowcell}_{refid}.bam"
    output:
        "samtools_stats_{refid}/{giab_id}_{flowcell}_{refid}.txt"
    log:
        "logs/samtools_stats/{giab_id}_{flowcell}_{refid}.log"
    threads: 1
    wrapper:
        "0.38.0/bio/samtools/stats"

rule bam_stats:
    input:
        "bams/{giab_id}_{flowcell}_{refid}.bam"
    output:
        "qc/{giab_id}_{flowcell}_{refid}.bam.stats.tsv.gz"
    log: "logs/bam_stats_{giab_id}_{flowcell}_{refid}.log"
    conda: "envs/bam_stats.yaml"
    shell:
        "python scripts/get_bam_stat.py {input} {output}"
    
rule multiqc:
    input:
        expand("samtools_stats_{{refid}}/{{giab_id}}_{flowcell}_{{refid}}.txt",
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


#rule get_hs37d5_phased_vars:
#    input: 
#        hs37d5 = FTP.remote(config["hs37d5"]["vars"]["{giab_id}"]),
#    output:
#        hs37d5 = "resources/hs37d5_{giab_id}.vcf.gz",
#        hs37d5_idx = "resources/hs37d5_{giab_id}.vcf.gz.tbi",
#    run:
#            ## Moving files to resources directory
#           shell("mv {input.hs37d5} {output.hs37d5}")
#            shell("tabix -p vcf {output.hs37d5}")

#rule get_GRCh38_phased_vars:
#    input: 
#        GRCh38 = FTP.remote(config["GRCh38"]["vars"]["{giab_id}"], keep_local = True)
#    output:
#        GRCh38 = "resources/GRCh38_{giab_id}.vcf.gz",
#        GRCh38_idx = "resources/GRCh38_{giab_id}.vcf.gz.tbi"
#    shell: """
#            ## Fixing sample names for GRCh38 to match bams
#	    ## This code is specifically for HG002 - many need to modify for other GIABIDs
#            zcat {input.GRCh38} \
#                | sed 's/26897$/HG002/' \
#                | awk '{{if($0~/^#/){{print $0}} else if(($7=="PASS") || ($7==".")){{print $0}}}}' \
#                | bgzip > {output.GRCh38}
#            tabix -p vcf {output.GRCh38}
#    """


rule phase_bams:
    input: 
        bam = "bams/{giab_id}_{flowcell}_{refid}.bam", 
        bamidx = "bams/{giab_id}_{flowcell}_{refid}.bam",
        ref = "resources/{refid}.fna",
        refidx = "resources/{refid}.fna.fai",
        var = "resources/{refid}_{giab_id}.vcf.gz",
        varidx = "resources/{refid}_{giab_id}.vcf.gz.tbi",
    output: "bams/{giab_id}_{flowcell}_{refid}_phased.bam"
    log: "logs/phase_bams_{giab_id}_{flowcell}_{refid}.log"
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
    output: "bams/{giab_id}_{flowcell}_phased_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"

## Combine bams --------------------------------------------------------------------
def get_bams(wildcards):
    gid_flowcells = flowcells.index[flowcells['giab_id'] == wildcards.giab_id]
    bams = expand("bams/{giab_id}_{flowcell}_{{refid}}_phased.bam", giab_id = wildcards.giab_id, flowcell = list(gid_flowcells))
    return(bams)

rule combine_bams:
    input: get_bams
    output: "{giab_id}_phased_{refid}.bam"
    log: "logs/combine_bams_{giab_id}_{refid}.log"
    params: "-f -O bam"
    threads: 8
    wrapper: "0.38.0/bio/samtools/merge"

rule index_combined:
    input: "{giab_id}_phased_{refid}.bam"
    output: "{giab_id}_phased_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"


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

## Running pipeline
rule all:
    input: 
        expand("bams/NA12878-minion-ul_{refid}.bam", refid = REFIDS), 
        expand("bams/NA12878-minion-ul_{refid}.bam.bai", refid = REFIDS),
        expand("qc/{flowcell}_{refid}.bam.stats.tsv.gz", 
               refid = REFIDS,
               flowcell = list(flowcells.index)),
        expand("phased/NA12878-minion-ul_{refid}.bam", refid = REFIDS), 
        expand("phased/NA12878-minion-ul_{refid}.bam.bai", refid = REFIDS)

        
def get_config(wildcards):
    return flowcells.loc[wildcards.flowcell,"config"]

## Aligning to reference -----------------------------------------------------------
def get_readgroup(wildcards):
    ## Format bam read group 
    run      = flowcells.loc[wildcards.flowcell, 'sample_id']
    model    = flowcells.loc[wildcards.flowcell, 'platform']
    lib      = flowcells.loc[wildcards.flowcell, 'sample_id']
    gpy_cfg  = flowcells.loc[wildcards.flowcell,'config']
    flowcell = flowcells.loc[wildcards.flowcell, 'type']
    kit      = flowcells.loc[wildcards.flowcell, 'kit']
    sample   = flowcells.loc[wildcards.flowcell, 'sample_id']
    
    read_group =   (
        f"@RG\\tID:{wildcards.flowcell}\\t"
        f"PL:nanopore\\t"
        f"PM:{model}\\t"
        f"LB:{lib}\\t"
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

rule map_reads:
    input:
        ref="resources/{refid}.fna", 
        refidx = "resources/{refid}.fna.fai",
        fastq="fastq/{flowcell}.fastq.gz"
    output: temp("sams/{flowcell}_{refid}.sam")
    params: read_group=get_readgroup, threads=8
    conda: "envs/map_reads.yaml"
    shell: """
        minimap2 -t {params.threads} --MD -aL -z 600,200 -x map-ont \
                -R \'{params.read_group}\' {input.ref} {input.fastq} > {output}
    """

rule sort_bam:
    input: 
        ref="resources/{refid}.fna", 
	sam="sams/{flowcell}_{refid}.sam"
    output: "bams/{flowcell}_{refid}.bam",
    params: mem=4, threads=8
    conda: "envs/map_reads.yaml"
    shell: """
        samtools sort \
		-m {params.mem}G -@{params.threads} \
                --reference {input.ref} \
		-o {output} \
		{input.sam}
    """

rule index_bam:
    input: "bams/{flowcell}_{refid}.bam"
    output: "bams/{flowcell}_{refid}.bam.bai"
    conda: "envs/map_reads.yaml"
    shell: " samtools index {output}"
    
## BAM QC 
rule bam_stats:
    input:
        "bams/{flowcell}_{refid}.bam"
    output:
        "qc/{flowcell}_{refid}.bam.stats.tsv.gz"
    conda: "envs/bam_stats.yaml"
    shell:
        "python scripts/get_bam_stat.py {input} {output}"


## Phasing reads --------------------------------------------------------------------

rule get_hs37d5_phased_vars:
    input: "resources/sp_v37.7.0.NA12878.vcf.gz"
    output:
        hs37d5 = "resources/hs37d5.vcf.gz",
        hs37d5_idx = "resources/hs37d5.vcf.gz.tbi",
    shell: """
            ## Fixing sample names for GRCh37 to match bams
            bcftools view  {input} \
                | sed 's/NA12878$/HG001/' \
                | bgzip > {output.hs37d5}
            tabix -p vcf {output.hs37d5}
    """

rule get_GRCh38_phased_vars:
    input: "resources/sp_v38.1.0.NA12878.vcf.gz"
    output:
        GRCh38 = "resources/GRCh38.vcf.gz",
        GRCh38_idx = "resources/GRCh38.vcf.gz.tbi"
    shell: """
            ## Fixing sample names for GRCh38 to match bams
            bcftools view {input} \
                | sed 's/NA12878$/HG001/' \
                | bgzip > {output.GRCh38}
            tabix -p vcf {output.GRCh38}
    """


rule phase_bams:
    input: 
        bam = "bams/NA12878-minion-ul_{refid}.bam", 
        bamidx = "bams/NA12878-minion-ul_{refid}.bam.bai",
        ref = "resources/{refid}.fna",
        refidx = "resources/{refid}.fna.fai",
        var = "resources/{refid}.vcf.gz",
        varidx = "resources/{refid}.vcf.gz.tbi",
    output: "phased/{flowcell}_{refid}.bam"
    params: whatshap = config["whatshap"]
    conda: "envs/map_reads.yaml"
    shell: """
         {params.whatshap} haplotag \
             -o {output} -r {input.ref} \
             {input.var} {input.bam}
     """

rule index_phased:
    input: rules.phase_bams.output
    output: "phased/{flowcell}_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"

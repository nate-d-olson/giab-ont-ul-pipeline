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

## Using signularity 
singularity: "docker://continuumio/miniconda3:4.4.10"

## Running pipeline
rule all:
    input: 
        expand("combined_{refid}.bam", refid = REFIDS), 
        expand("combined_{refid}.bam.bai", refid = REFIDS),
        "combined.fastq.gz",
	"combined.sequencing_summary.txt.gz",
        expand("combined_{refid}_qc.txt", refid = REFIDS),
        expand("qc/{flowcell}_{refid}.bam.stats.tsv.gz", 
               refid = REFIDS,
               flowcell = list(flowcells.index)),
         expand("phased_{refid}.bam", refid = REFIDS), 
         expand("phased_{refid}.bam.bai", refid = REFIDS)


def get_config(wildcards):
    return flowcells.loc[wildcards.flowcell,"config"]

## Basecalling --------------------------------------------------------------------
rule basecalling:
    input: "../../raw/{flowcell}"
    output: directory("fastq/{flowcell}")
    params: get_config
    shell: 'sw/ont-guppy/bin/guppy_basecaller -x cuda:0 -r \
                --num_callers 8 \
                --gpu_runners_per_device 4 \
                --chunks_per_runner 1664 \
                --input_path {input} \
                --save_path {output} \
                --records_per_fastq 0 \
                --config {params}'

rule combine_and_compress_fastq:
    input: rules.basecalling.output
    output: "fastq/{flowcell}.fastq.gz"
    shell: 'cat {input}/*fastq | gzip > {output}'

## Combining fastq and sequencing summaries ---------------------------------------

rule combine_fastq:
    input: expand("fastq/{flowcell}.fastq.gz", flowcell = list(flowcells.index))
    output: "combined.fastq.gz"
    shell: "cat {input} > {output}"

rule combine_seq_summary:
    input: expand("fastq/{flowcell}/sequencing_summary.txt", flowcell = list(flowcells.index))
    output: "combined.sequencing_summary.txt.gz"
    run: 
        combined_tsv = pd.concat([pd.read_csv(f , delimiter='\t') for f in input])
        combined_tsv.to_csv( output[0], sep = "\t", 
                             index = False, compression = "gzip")

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
    
    read_group =   (
        f"@RG\\tID:{wildcards.flowcell}\\t"
        f"PU:{wildcards.flowcell}\\t"
        f"PL:nanopore\\t"
        f"PM:{model}\\t"
        f"LB:{lib}\\t"
        f"DT:{date}\\t"
        f"PG:guppy-{config['guppy_version']}-{gpy_cfg}\\t"
        f"DS:Flowcell={flowcell},kit={kit}\\t"
        f"SM:HG002"
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
        fastq=rules.combine_and_compress_fastq.output
    output: 
        bam="bams/{flowcell}_{refid}.bam",
        bai="bams/{flowcell}_{refid}.bam.bai"
    params: read_group=get_readgroup, mem=4, threads=4, temp_id="{flowcell}_{refid}"
    conda: "envs/map_reads.yaml"
    shell: """    
        minimap2 -t {params.threads} -a -L -z 600,200 -x map-ont \
                -R \'{params.read_group}\' {input.ref} {input.fastq} > sams/{params.temp_id}.sam
        
        samtools sort \
		-m {params.mem}G -@{params.threads} \
                --reference {input.ref} \
		-o {output.bam} \
		sams/{params.temp_id}.sam

        samtools index {output.bam}
    """

#    shell: """
#        minimap2 -t {params.threads} -a -L -z 600,200 -x map-ont \
#                -R \'{params.read_group}\' {input.ref} {input.fastq} \
#            | samtools sort --verbosity 5 -m {params.mem}G -@{params.threads} \
#                -O bam --reference {input.ref} > {output}
#        samtools index {output}
#    """
    
## Flowcell BAM QC ----------------------------------------------------------------
rule samtools_stats:
    input: "bams/{flowcell}_{refid}.bam"
    output: "samtools_stats_{refid}/{flowcell}_{refid}.txt"
    log: "logs/samtools_stats/{flowcell}_{refid}.log"
    wrapper: "0.38.0/bio/samtools/stats"

rule bam_stats:
    input: "bams/{flowcell}_{refid}.bam"
    output: "qc/{flowcell}_{refid}.bam.stats.tsv.gz"
    conda: "envs/bam_stats.yaml"
    shell: "python scripts/get_bam_stat.py {input} {output}"

## Combine bams --------------------------------------------------------------------
rule combine_bams:
    input: expand("bams/{flowcell}_{{refid}}.bam", flowcell = list(flowcells.index))
    output: "combined_{refid}.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule index_combined:
    input: "combined_{refid}.bam"
    output: "combined_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"


## Combined BAM QC
rule combined_bam_stats:
    input: "combined_{refid}.bam"
    output: 
        txt="combined_{refid}_qc.txt"
    conda: "envs/bam_stats.yaml"
    shell: "python scripts/quick_qc.py {input} > {output.txt}"

## Phasing reads --------------------------------------------------------------------
rule get_hs37d5_phased_vars:
    input: 
        hs37d5 = FTP.remote(config["hs37d5"]["vars"]),
    output:
        hs37d5 = "resources/hs37d5.vcf.gz",
        hs37d5_idx = "resources/hs37d5.vcf.gz.tbi",
    run:
            ## Moving files to resources directory
            shell("mv {input.hs37d5} {output.hs37d5}")
            shell("tabix -p vcf {output.hs37d5}")

rule get_GRCh38_phased_vars:
    input: 
        GRCh38 = FTP.remote(config["GRCh38"]["vars"], keep_local = True)
    output:
        GRCh38 = "resources/GRCh38.vcf.gz",
        GRCh38_idx = "resources/GRCh38.vcf.gz.tbi"
    shell: """
            ## Fixing sample names for GRCh38 to match bams
            zcat {input.GRCh38} \
                | sed 's/26897$/HG002/' \
                | awk '{{if($0~/^#/){{print $0}} else if(($7=="PASS") || ($7==".")){{print $0}}}}' \
                | bgzip > {output.GRCh38}
            tabix -p vcf {output.GRCh38}
    """

rule phase_bams:
    input: 
        bam = "bams/{flowcell}_{refid}.bam", 
        bamidx = "bams/{flowcell}_{refid}.bam.bai",
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

## Combine Phased Reads -------------------------------------------------------------
rule combine_phased:
    input: expand("phased/{flowcell}_{{refid}}.bam", flowcell = list(flowcells.index))
    output: "phased_{refid}.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule index_combined_phased:
    input: rules.combine_phased.output
    output: "phased_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"

## Combined BAM QC
rule phased_bam_stats:
    input: "phased_{refid}.bam"
    output:"phased_{refid}_qc.txt"
    conda: "envs/bam_stats.yaml"
    shell: "python scripts/quick_qc.py {input} > {output}"

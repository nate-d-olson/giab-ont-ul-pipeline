## Pipeline dependencies - samtools, bcftools, guppy, and whatshap
import pandas as pd

## For ftp files
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

## Pipeline Config
configfile: "config.yaml"

## Reference Genomes
REFIDS=['hs37d5','GRCh38']

## Variables
outdir=config["giab_id"]

## Loading flowcell metadata
flowcells = pd.read_csv("flowcell_metadata.csv").set_index("flowcell_id")

## Using signularity
singularity: "docker://continuumio/miniconda3:4.4.10"


## Running pipeline
rule all:
    input:
        expand(outdir + "qc/{flowcell}_{refid}.bam.stats.tsv.gz",
               refid = REFIDS,
               flowcell = list(flowcells.index)),
        expand(outdir + "phased_{refid}.bam", refid = REFIDS),
        expand(outdir + "phased_{refid}.bam.bai", refid = REFIDS),
        expand(outdir + "phased_{refid}_qc.txt", refid = REFIDS),
        expand(outdir + "phased_{refid}.bam.stats.tsv.gz", refid = REFIDS),
        expand(outdir + "phased_{refid}.flagstat.txt", refid = REFIDS),
        expand(outdir + "phased_{refid}.stats.txt", refid = REFIDS)


def get_config(wildcards):
    return flowcells.loc[wildcards.flowcell,"config"]

## Combining fastq and sequencing summaries ------------------------------------
rule combine_fastq:
    input: expand("fastq/{flowcell}.fastq.gz", flowcell = list(flowcells.index))
    output: outdir + "combined.fastq.gz"
    shell: "cat {input} > {output}"

rule qc_combined_fastq:
    """
    Input is the reads in directory output is info about reads
    """
    input: outdir + "combined.fastq.gz"
    output: outdir + "combined.read_stats.txt"
    message: "Calculating read coverage statitics for: {input}",
    params:
        read_stat_script = "scripts/rawcoverage.py",
    threads: 12
    conda: "envs/rawcoverage.yaml"
    shell:
        """
        python {params.read_stat_script} -i {input} -o {output} -t {threads}
        """

rule combine_seq_summary:
    input: expand(outdir + "fastq/{flowcell}/sequencing_summary.txt", flowcell = list(flowcells.index))
    output: outdir + "combined.sequencing_summary.txt.gz"
    run:
        combined_tsv = pd.concat([pd.read_csv(f , delimiter='\t') for f in input])
        combined_tsv.to_csv( output[0], sep = "\t",
                             index = False, compression = "gzip")

## Aligning to reference -------------------------------------------------------
def get_readgroup(wildcards):
    ## Format bam read group
    run      = flowcells.loc[wildcards.flowcell, 'sample_id']
    model    = flowcells.loc[wildcards.flowcell, 'platform']
    lib      = flowcells.loc[wildcards.flowcell, 'sample_id']
    date     = flowcells.loc[wildcards.flowcell, 'date']
    gpy_cfg  = flowcells.loc[wildcards.flowcell,'config']
    flowcell = flowcells.loc[wildcards.flowcell, 'type']
    kit      = flowcells.loc[wildcards.flowcell, 'kit']
    sample   = config["giab_id"]

    read_group =   (
        f"@RG\\tID:{wildcards.flowcell}\\t"
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

rule map_reads:
    input:
        ref="resources/{refid}.fna",
        refidx = "resources/{refid}.fna.fai",
        fastq=rules.combine_and_compress_fastq.output
    output:
        bam= outdir + "bams/{flowcell}_{refid}.bam",
        bai= outdir + "bams/{flowcell}_{refid}.bam.bai"
    params: read_group=get_readgroup, mem=12, threads=4
    conda: "envs/map_reads.yaml"
    shell: """
        minimap2 -t {params.threads} -a -L -z 600,200 -x map-ont \
                -R \'{params.read_group}\' {input.ref} {input.fastq} \
            | samtools sort -m {params.mem}G -@{params.threads} \
                -O bam --reference {input.ref} > {output.bam}
        samtools index {output.bam}
    """

## Flowcell BAM QC -------------------------------------------------------------
rule samtools_stats:
    input: outdir + "bams/{flowcell}_{refid}.bam"
    output: outdir + "samtools_stats_{refid}/{flowcell}_{refid}.txt"
    log: "logs/samtools_stats/{flowcell}_{refid}.log"
    wrapper: "0.38.0/bio/samtools/stats"

rule bam_stats:
    input: "bams/{flowcell}_{refid}.bam"
    output: outdir + "qc/{flowcell}_{refid}.bam.stats.tsv.gz"
    conda: "envs/bam_stats.yaml"
    shell: "python scripts/get_bam_stat.py {input} {output}"

## Combine bams ----------------------------------------------------------------
rule combine_bams:
    input: expand(outdir + "bams/{flowcell}_{{refid}}.bam",
    flowcell = list(flowcells.index))
    output: outdir + "combined_{refid}.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule index_combined:
    input: "combined_{refid}.bam"
    output: outdir + "combined_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"


## Phasing reads ---------------------------------------------------------------
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
        GRCh38 = "resources/" + config["giab_id"] + "GRCh38.vcf.gz",
        GRCh38_idx = "resources/" + config["giab_id"] + "GRCh38.vcf.gz.tbi"
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
        bam = outdir + "bams/{flowcell}_{refid}.bam",
        bamidx = outdir + "bams/{flowcell}_{refid}.bam.bai",
        ref = "resources/{refid}.fna",
        refidx = "resources/{refid}.fna.fai",
        var = "resources/" + config["giab_id"] + "{refid}.vcf.gz",
        varidx = "resources/" + config["giab_id"] + "{refid}.vcf.gz.tbi",
    output: outdir + "phased/{flowcell}_{refid}.bam"
    params: whatshap = config["whatshap"]
    conda: "envs/map_reads.yaml"
    shell: """
        {params.whatshap} haplotag \
            -o {output} -r {input.ref} \
            {input.var} {input.bam}
    """

rule index_phased:
    input: rules.phase_bams.output
    output: outdir + "phased/{flowcell}_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"

## Combine Phased Reads --------------------------------------------------------
rule combine_phased:
    input: expand("phased/{flowcell}_{{refid}}.bam",
                flowcell = list(flowcells.index))
    output: outdir + "phased_{refid}.bam"
    params: "-f -O bam"
    threads: 3
    wrapper: "0.38.0/bio/samtools/merge"

rule index_combined_phased:
    input: rules.combine_phased.output
    output: outdir + "phased_{refid}.bam.bai"
    wrapper: "0.38.0/bio/samtools/index"

## Combined BAM QC
rule phased_bam_stats:
    input: outdir + "phased_{refid}.bam"
    output: outdir + "phased_{refid}_qc.txt"
    conda: "envs/bam_stats.yaml"
    shell: "python scripts/quick_qc.py {input} > {output}"

rule phased_bam_full_stats:
    input: outdir + "phased_{refid}.bam"
    output: outdir + "phased_{refid}.bam.stats.tsv.gz"
    conda: "envs/bam_stats.yaml"
    shell: "python scripts/get_bam_stat.py {input} {output}"


rule phased_flagstat:
    input: outdir + "phased_{refid}.bam"
    output: outdir + "phased_{refid}.flagstat.txt"
    wrapper: "0.49.0/bio/samtools/flagstat"


rule phased_stats:
    input: outdir + "phased_{refid}.bam"
    output: outdir + "phased_{refid}.stats.txt"
    wrapper: "0.49.0/bio/samtools/stats"


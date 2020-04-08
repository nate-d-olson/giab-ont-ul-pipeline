#!/usr/bin/bash -l
#SBATCH -n 4  # number of cores
#SBATCH --mem 8g # memory pool for all cores
#SBATCH --job-name=guppy-lo
#SBATCH --time=2-00:00:00
#SBATCH -p gpu,owners
#SBATCH --gres gpu:2
#SBATCH -C GPU_GEN:VLT
#SBATCH --mail-type=ALL


## Loading required modules
module load python/3.6.1 R

FLOWCELLID=PAC17886
GUPPY=sw/ont-guppy/bin/guppy_basecaller
F5DIR=../../raw_compressed/${FLOWCELLID}
FQDIR=fastq/${FLOWCELLID}

## Get GPU information
nvidia-smi

## Running guppy
${GUPPY} -x "cuda:all" -r \
    --num_callers 8 \
    --gpu_runners_per_device 4 \
    --chunks_per_runner 1664 \
    --input_path ${F5DIR} \
    --save_path ${FQDIR} \
    --records_per_fastq 0 \
    --config dna_r9.4.1_450bps_hac_prom.cfg
#    --resume \

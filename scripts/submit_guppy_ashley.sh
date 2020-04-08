#!/usr/bin/bash -l
#SBATCH -n 4  # number of cores
#SBATCH --mem 8g # memory pool for all cores
#SBATCH --array=0
#SBATCH --job-name=guppy-ash
#SBATCH --output=array_%A_%a.out
#SBATCH --error=array_%A_%a.err
#SBATCH --time=8:00:00
#SBATCH -p gpu,owners
#SBATCH --gres gpu:2
#SBATCH -C GPU_GEN:VLT
#SBATCH --mail-type=ALL


## Loading required modules
module load python/3.6.1 R

FLOWCELLID=PAD11290
GUPPY=sw/ont-guppy/bin/guppy_basecaller
F5DIR=../../raw_compressed/${FLOWCELLID}/${FLOWCELLID}
FQDIR=fastq/${FLOWCELLID}

## Get GPU information
nvidia-smi

## Index id
i=$SLURM_ARRAY_TASK_ID

## Running guppy
# for i in {1..10}; do    
    ${GUPPY} -x "cuda:all" -r \
        --num_callers 8 \
        --gpu_runners_per_device 4 \
        --chunks_per_runner 1664 \
        --input_path ${F5DIR}_${i} \
        --save_path ${FQDIR}_${i} \
        --records_per_fastq 0 \
        --config dna_r9.4.1_450bps_hac_prom.cfg
#done

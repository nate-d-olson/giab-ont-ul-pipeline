# Pipeline for processing JIMB UL ONT dataset 

Can use helper script `pipe.sub` for pipeline job submission.

Start by activating miniconda
```
source ~/Miniconda/bin/activate
```

Create conda environment with Snakemake and pipeline global dependencies defined in `environment.yaml`.

```
conda env create -n ont-ul-pipe
```

After creating the pipeline environment, activate the environment using `conda activate ont-ul-pipe`

Run pipeline using a command similar to the one below, code after `--cluster-config` only for running on a cluster, `pipe.sub` is a helper script for pipeline job submission.
```
snakemake --use-conda -j 999 --cluster-config cluster.json \
	--cluster "sbatch -J {cluster.name} \
			-p {cluster.partition} \ 
			-n {cluster.n}  \
			--mem {cluster.mem} \
			-t {cluster.time} \
			{cluster.resources} \
			-e {cluster.error} \
			-o {cluster.output} \
			--mail-type ALL \
			--profile ALL"
```

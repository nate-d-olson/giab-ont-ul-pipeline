# Pipeline for processing JIMB UL ONT dataset 

Can use helper script `pipe.sub` for pipeline job submission.

start by activating miniconda with snakemake installed
`source ~/Miniconda/bin/activate`

run pipeline using 
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

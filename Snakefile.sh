#!/bin/sh

module load python3/3.10.2 slurm

DATE=$(date +%y%m%d)
mkdir -p logs_${DATE}

snakemake --cores=2 --unlock
sbcmd="sbatch --partition=cgrq --mem 64g --cpus-per-task={threads} --output=logs_${DATE}/snakejob_%j.out"
snakemake --use-conda -pr --cluster "$sbcmd" --keep-going --rerun-incomplete --jobs 30 --latency-wait 60

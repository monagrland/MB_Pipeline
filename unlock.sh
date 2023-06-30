#!/bin/bash
# A shell script to unlock the directories
eval "$(conda shell.bash hook)"
conda env create -n mb_snakemake --file "$PWD/envs/mb_snakemake.yaml" #creates the conda environment if you have the correct permissions and if it doesnt exist yet.
THREADS=($(grep "threads:" $1 | sed s/"threads: "/""/))
conda activate mb_snakemake
echo "Activated conda environment."
snakemake --use-conda --cores $THREADS --conda-frontend conda --configfile $1 --unlock
echo "Finished unlocking."

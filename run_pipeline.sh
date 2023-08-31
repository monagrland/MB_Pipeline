#!/bin/bash
# A shell script to start the pipeline with the specified config file in the correct conda environment
eval "$(conda shell.bash hook)"
conda env create -n mb_snakemake --file "$PWD/envs/mb_snakemake.yaml" #creates the conda environment if you have the correct permissions and if it doesnt exist yet.
STARTTIME=$(date +%s)
THREADS=($(grep "threads:" $1 | sed s/"threads: "/""/))
conda activate mb_snakemake
echo "Activated conda environment."
echo "Starting the pipeline using $THREADS threads."
snakemake --use-conda --cores $THREADS --conda-frontend conda --notemp --configfile $1
ENDTIME=$(date +%s)
echo "It took $(($ENDTIME - $STARTTIME)) seconds to finish this run."

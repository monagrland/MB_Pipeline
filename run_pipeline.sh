#!/bin/bash
# A shell script to start the pipeline with the specified config file in the correct conda environment
eval "$(conda shell.bash hook)"
conda env create -n mb_snakemake --file "$PWD/envs/mb_snakemake.yaml" #creates the conda environment if you have the correct permissions and if it doesnt exist yet.
STARTTIME=$(date +%s)
conda activate mb_snakemake
echo "activated conda environment"
snakemake --use-conda --cores 8 --conda-frontend conda --configfile $1
ENDTIME=$(date +%s)
echo "It took $(($ENDTIME - $STARTTIME)) seconds for this Pipeline"

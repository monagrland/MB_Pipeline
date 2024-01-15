#!/bin/bash

# Assumes that conda is in path, and that mamba and conda are in base environment
eval "$(conda shell.bash hook)"
# Test for snakemake environment in current folder
mamba list -p ./snakemake_8 &> /dev/null
RETVAL=$?
if [[ $RETVAL == 1 ]]
then
  echo "Creating environment with Snakemake at ./snakemake_8"
  mamba env create -p ./snakemake_8 --file "$PWD/workflow/envs/snakemake.yaml"
else
  echo "Environment with Snakemake found at ./snakemake_8"
fi
STARTTIME=$(date +%s)
CORES=12
conda activate ./snakemake_8
echo "Activated conda environment."
echo "Starting the pipeline using $THREADS threads."
snakemake --sdm conda --cores $CORES --conda-frontend mamba --notemp --configfile config/config.yaml
ENDTIME=$(date +%s)
echo "It took $(($ENDTIME - $STARTTIME)) seconds to finish this run."

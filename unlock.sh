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
conda activate ./snakemake_8
echo "Activated conda environment."
snakemake --sdm conda --conda-frontend mamba --notemp --configfile config/config.yaml --unlock $@
echo "Unlocked"

#!/bin/bash

source activate snakemake # activate your environment containing snakemake
pipdir=/path/to/scNanoHi-C/pipeline
njob=100

read -p "Run scNanoHi-C snakemake pipeline? [y/dry/n]:" option

if [ $option == "y" ];then
    snakemake --profile slurm --jobs $njob \
    -s $pipdir/filter_contacts.smk \
    --configfile ./config.yaml -k
elif [ $option == "dry" ];then
    snakemake --profile slurm --jobs $njob \
    -s $pipdir/filter_contacts.smk \
    --configfile ./config.yaml -k -n -p
else
    echo "Canceled!"
    exit 0
fi
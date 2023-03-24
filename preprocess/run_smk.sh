#!/bin/bash
base_dir="/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen"
source ${base_dir}/source/ljs_bashrc.sh
source activate snakemake
pipdir=${base_dir}/Github/Pipeline/scHic
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
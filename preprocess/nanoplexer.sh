#!/bin/bash
#SBATCH -p cn-long
#SBATCH -c 16
#SBATCH -J LJS_nanoplx
#SBATCH -o nanoplx.%j.out
#SBATCH -e nanoplx.%j.err
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchoucnl
source activate nanopore
nanoplexer \
	-b barcode.fa \
	-t 16 \
	-p raw_data \
	pass.fastq.gz

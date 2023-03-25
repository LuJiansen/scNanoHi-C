#!/bin/bash
#SBATCH -p cn-long
#SBATCH -c 20
#SBATCH -J LJS_QUAST
#SBATCH -o quast.%j.out
#SBATCH -e quast.%j.err
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchoucnl

source activate HiC
asm=scaffolds_FINAL.fasta
ref=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/database/GRCh38_ref/GRCh38.fa
gtf=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/database/hg38_ensembl_104/Homo_sapiens.GRCh38.104.gtf
prefix=quast_flt

quast $asm \
	-r $ref \
	-g $gtf \
	-o $prefix \
	-x 5000 \
	-m 10000 \
	-i 10000

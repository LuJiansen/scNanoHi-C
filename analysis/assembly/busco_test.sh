#!/bin/bash
#SBATCH -p cn-long
#SBATCH -c 20
#SBATCH -N 1
#SBATCH -J LJS_BUSCO
#SBATCH -o busco.%j.out
#SBATCH -e busco.%j.err
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchoucnl

source activate HiC
download_path=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/database/busco_downloads
lineage=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/database/busco_downloads/lineages/vertebrata_odb10
busco -i scaffolds_FINAL.fasta -l $lineage \
	-o busco -m genome --offline \
	--download_path $download_path

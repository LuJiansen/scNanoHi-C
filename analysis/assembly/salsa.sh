#!/bin/bash
#SBATCH -p cn-long
#SBATCH -c 10
#SBATCH -J LJS_SALSA
#SBATCH -o salsa.%j.out
#SBATCH -e salsa.%j.err
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchoucnl

fa=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/project/scHic/assembly/hg002_10_nano_flye/assembly.fasta
fai=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/project/scHic/assembly/hg002_10_nano_flye/assembly.fasta.fai
gfa=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/project/scHic/assembly/hg002_10_nano_flye/assembly_graph.gfa
enzyme=GATC #MboI

bed=$1
prefix=$2

echo "fasta: $fa"
echo "index: $fai"
echo "graph: $gfa"
echo "enzyme site: $enzyme"
echo "alignments: $bed"
echo "output dir: $prefix"

source activate LJS_py2

run_pipeline.py -a $fa \
	-l $fai \
	-b $bed \
	-e $enzyme \
	-o scaffolds_${prefix} \
	-m yes \
	-g $gfa

# scNanoHi-C
  scNanoHi-C is a single-cell long-read concatemer sequencing method which could be used to investigate the higher-order 3D genomic structures in individual cells. This repository contains snakemeke workflows to preprocess data from scNanoHi-C

# Usage
1. data_processing.smk is used to demultiplex and remove sequencing adaptors
2. filter_contacts.smk is used to remove artifacts of single-cell Hi-C, perform basic quality control and reconstruct single-cell 3D genomic structure.

# README

Last edited time: March 24, 2023 3:47 PM

## Introduction

scNanoHi-C is a single-cell long-read concatemer sequencing method which could be used to investigate the higher-order 3D genomic structures in individual cells. This repository contains snakemeke workflows to preprocess scNanoHi-C data and other codes for analyses in our paper.

## Usage

### Demultiplex and remove adapter sequences

- Input files:
    - `pass.fastq.gz` the output sequencing data from Nanopore platform;
    - `bc_index` a tab-separated file containing PCR barcodes (1-96) in the first column and an optional library name in second column as follows:
        
        ```
        $ head bc_index
        65
        66
        67
        68
        ```
        
- Output files:
    - `raw_data/`  containing raw fastq files for each single cell
    - `trim/` containing trimmed fastq files for each single cell, can be used as the input files for [Pore-C-Snakemake](https://github.com/nanoporetech/Pore-C-Snakemake)
    - `raw_fastq.stats` statistics of raw fastq files
    - `trim_fastq.stats` statistics of trimmed fastq files
- usage:
    
    put `pass.fastq.gz` and  `bc_index` together with `data_processing.smk` into the working directory and run snakemake via:
    
    ```bash
    snakemake -s data_processing.smk -j 10
    ```
    

### Run Pore-C-Snakemake to generate concatemers

trimmed fastq files for single cell were used to run default Pore-C-Snakemake workflow, see details in [Pore-C-Snakemake](https://github.com/nanoporetech/Pore-C-Snakemake).

### Run single-cell filtration and quality control pipeline

The main script of this step is `filter_contacts.smk` , which could be  run with the shell script `run_smk.sh` . This step is consist of following sections

- **remove artifact contacts in single cells**
    
    To remove artifact contacts, concatemers from each single cell were firstly decomposed into VPCs. Five types of artifacts were removed sequentially:
    
    1. Adjacent contacts: contacts assigned to adjacent restriction fragments;
    2. Close contacts: contacts with two alignments separating from less than 1000 bp in genomic distance;
    3. Duplicates: we designated contacts from the same pair of restriction fragments as duplicates. In order to remove PCR duplicates and preserve the high-order structure of concatemers as much as possible, only one contact from the concatemer with the highest cardinality in each duplicates set was selected to remain.
    4. Promiscuous contacts: contacts from restriction fragments which involved in more than 10 interactions.
    5. Isolated contacts: contacts which no other contact within 1Mb distance.
    
    More information about the evaluation of these filters could be found in the Methods of our paper.
    
- **reconstruct single-cell 3D models**
    
    The haplotype-tagged virtual pair-wise contacts from scNanoHi-C were used for haplotype imputation and single-cell 3D genome structure model reconstruction by [hickit](https://github.com/lh3/hickit) and [dip-C](https://github.com/tanlongzhi/dip-c) packages.The high-depth data (24 cells per run) of scNanoHi-C were recommended for 3D genome construction.
    
- **generate quality control metric for downstream analysis**
- **calculated the single-cell A/B compartment values (scA/B)**
    
    The scA/B values were calculated through the dip-C package in both 2D and 3D mode, which could be used for dimensionality reduction and cell clustering. The recommended resolution is 1 Mb.
    
- Input files:
    - `run_smk.sh` : the main script.
    - results from Pore-C pipeline: including the `.parquet` files and concatemer summary file of each cell.
    - `config.yaml` : the configure files containing the parameters and location of files used in the pipeline.
    - `sample` : a file containing the name of cells in each row.
- main output files:
    - `*.contacts.clean.parquet` : the files containing filtered contacts for each cell.
    - `*_stats_summary.csv`: files containing the quality control information of all cells.
    - `*.n.cif` : results of single-cell 3D models.
    - `RMSD_summary.txt` : summary of RMSD of all models.
    - `cpg_2d_*_b1m_summary.csv` : matrix of 2D scA/B values in all cells.
    
    These outputs could be optional by adjusting the `rule all` in `filter_contacts.smk` .
    
- usage:
    
    put the `run_smk.sh` , `config.yaml` and `sample` files into the `merged_contacts` directory of Pore-C output, and modify these files accordingly then run as following:
    
    ```bash
    $ sh run_smk.sh
    Run scNanoHi-C snakemake pipeline? [y/dry/n]:
    # type 'y' to run the pipeline directly
    # type 'dry' to display what would be done
    # type 'n' to cancel the program
    ```
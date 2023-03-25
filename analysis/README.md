# Analysis of scNanoHi-C

Last edited time: March 25, 2023 5:18 PM

This folder contains the scripts used for analyses in scNanoHi-C, which including following sections:

### assembly

This section containing scripts for scaffolding single-cell *de novo* assembly with scNanoHi-C.  We used the [SALSA2](https://github.com/marbl/SALSA) pipeline to performed assembly scaffolding. In detail, the scNanoHi-C data were first mapped to the single-cell draft assemblies using scNanoHi-C snakemake pipeline, and the output filtered contacts of each cells were converted and merged into salsa bed format. Then the scaffolding was performed using the `run_pipeline.py` script of SALSA2.

- `get_salsa_bed.smk` : the Snakemake file to converted and merged scNanoHi-C output into SALSA2 format.
- `run_salsa.sh` : the script to run SALSA2 pipeline
- `run_quast.sh` : calculate QUAST metrics
- `run_busco.sh` : calculate BUSCO metrics

### clustering

- `scAB_clustering.R`: script for the single-cell dimensionality reduction and clustering using the single-cell A/B compartment (scA/B) values.
- `DCWs_enrichments.R` : script for identifying the differential conmpartmentalized windows (DCWs) between cell populations and performing the enrichment analyses of DCWs to the annotations of given cell types.

### detection_score

Scripts for calculating the single-cell detection scores of genomic structures based on the method descripted in the [scSPRITE](https://doi.org/10.1038/s41587-021-00998-1) paper.

- `Territories_detection_score.R`: script for calculated the territories detection scores among each pair of chromosomes.
- `Compartment_detection_score.R` : script for calculating the compartment detection scores of each A/B compartment switch triplets, including A-B-A and B-A-B.
- `Territories_detection_score.R`: script for calculating the TAD detection scores in up and down-stream 1 Mb of all strong TAD boundaries identified in merged scNanoHi-C data.

### ecDNA

Scripts for calculating the copy-number adjusted normalized trans-chromosomal interaction frequency (adjnTIF) with scNanoHi-C in cells carrying ecDNAs, using the similar methods described in [previous works](https://doi.org/10.1038/s41422-020-00466-6).

- `cal_adjnTIF.smk`: the Snakemake file to calculate the adjnTIF for each cells based on the mcool files of scNanoHi-C contacts.
- `cal_adjnTIF.R`: R script for calculating adjnTIF values.
- `cool2matrix`: script for converting the mcool files into 3-column tab-separated files each row of which represented a contact with the format `bin1_id bin2_id count` . For example:
    
    ```bash
    $ zcat example_mtx.gz | head
    0       106   1
    0       151   5
    0       153   1
    1       180   2
    1       234   4
    1       235   1
    ```
    

### genomic_variations

The scripts for identifying the CNV and SV in scNanoHi-C data. In brief, [NeoloopFinder](https://github.com/XiaoTaoWang/NeoLoopFinder) and [hic_breakerfinder](https://github.com/dixonlab/hic_breakfinder) were used to identify CNV and SV, respectively.

- `NeoloopFinder_CNV.smk`: the Snakemake file of the CNV calling pipeline for scNanoHi-C.
- `chronly.R` & `frags_flt.R`: scripts used in the `NeoloopFinder_CNV.smk` to filter non-chromosomal contigs in the contacts files and the restriction fragments.
- `hic_breakfinder_SV.smk`: the Snakemake file of hic_breakfinder SV calling pipeline for scNanoHi-C.
- `salsaBed2bam.sh`: script used in the `hic_breakfinder_SV.smk` to convert the input bed into bam used in hic_breakfinder.

### high-order_ratio

Scripts for the analyses about the ratio of high-order cnocatemers and cells in scNanoHi-C. The high-order ratios were defined to assess the probability of a genomic region being involved in multiway interactions in a given cell population.

- `high_order_ratios.R` : script for calculating the high-order ratio of concatemers and cells in scNanoHi-C.
- `enrichment_of_monomers.R` : enrichment analyses of monomers derived from high-order concatemers to the annotations of given cell types.
- `high_order_in_open_chromatin.R`: evalution of the relationships between high-order ratio and other factors.

### multiple_interaction

`multiple_regulatory_structures.R` : Script for identifying multiway regulatory structures supported by concatemers in scNanoHi-C, including the one promoter regulated by multiple enhancers and one enhancer simultaneously regulting multiple promoters. Enhancer-promoter relationships predicted by [ABC models](https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction) were used as reference.

### synergy

Scripts for identifying single-cell ‘synergies’ in scNanoHi-C. To assess the cooperativity of multiple genomic regions, we adapted the [Synergy and Chromunity](https://github.com/mskilab/chromunity) algorithms developed for Pore-C to identify synergies in our scNanoHi-C data, which includes candidate bin-sets identification and cooperativity test by Synergy model.

- `prepare_synergy_input.smk`: the Snakemake file containing the pipeline to prepare the input files for single-cell synergy analysis.
- `pq2monomer.R`:  script used in the `prepare_synergy_input.smk` to convert the parquet format of contacts into plain text file of monomers, for axample:
    
    ```bash
    $ zcat merged_monomer.txt.gz | head
    chrom   start   end     read_name       read_length     haplotype       read_idx        cell
    chr2    58089577        58090044        000354cf-3477-4941-ab59-872c9113ffbc    3837    2       3    GM_24.P94B10
    chr13   26678355        26679978        0004259b-4e1b-49eb-bcca-a3618ca34cd3    4520    1       5    GM_24.P94B10
    chr6    25794519        25795773        000bf675-1e9a-492f-9f4b-cdc7e08ef838    2494    1       17   GM_24.P94B10
    chr10   123686282       123686905       000c11cd-0744-4693-8952-22b63a00d51a    3012    -1      18   GM_24.P94B10
    chr10   120588390       120588650       000c11cd-0744-4693-8952-22b63a00d51a    3012    -1      18   GM_24.P94B10
    chr1    147275426       147275822       000dc096-39ba-4782-b472-1c5960d8b4c0    1892    -1      23   GM_24.P94B10
    chr1    147290186       147290760       000dc096-39ba-4782-b472-1c5960d8b4c0    1892    -1      23   GM_24.P94B10
    chr5    138394486       138395267       0011f175-9e2b-4785-a321-50450a25e7fe    2018    -1      31   GM_24.P94B10
    chr5    138533267       138533640       0011f175-9e2b-4785-a321-50450a25e7fe    2018    -1      31   GM_24.P94B10
    ```
    
- `synergy_funcs.R`: script containing the functions used in the single-cell synergy analysis.
- `run_all_synergy.R`: script to run single-cell synergy analysis around all 20 kb bins in each 2Mb windows.
- `run_RE_synergy.R`: script to run single-cell synergy analysis with all promoters and enhancers identified from the chromHMM annotations.
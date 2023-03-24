import os
import sys

sample_list = 'sample'
samples=[]
with open(sample_list,'r') as f:
	for line in f:
		samples.append(line.strip('\n'))
print(samples)
phases = config["phase"]
reference = config["reference"]
reps = config['reps']
resolutions = config['resolution']

chromsizes = '../refgenome/' + reference + '.rg.chromsizes'
fa = '../refgenome/' + reference + '.rg.fa'

wildcard_constraints:
    sample="[0-9A-Za-z]+"

print(reference)
print(phases)
print(reps)
print(resolutions)

rule all:
    input:
        expand("salsa/MboI_{sample}_{ref}_{phase}_filter.salsa.bed",
        sample = samples, ref = reference, phase = phases),
        expand("salsa/MboI_{sample}_{ref}_{phase}_filter.salsa.bam",
        sample = samples, ref = reference, phase = phases),
        "merged_salsa.bam",
        "merged_salsa.breaks.txt",

rule to_salsa_bed:
    input:
        contacts = "MboI_{sample}_{ref}_{phase}.contacts.clean.parquet",
    output:
        "salsa/MboI_{sample}_{ref}_{phase}_filter.salsa.bed",
    params:
        prefix = "salsa/MboI_{sample}_{ref}_{phase}_filter",
    log:
        "logs/{sample}_{ref}_{phase}.salsa.log",
    threads: 5,
    shell: """
        set +u; source activate HiC; set -u
        pore_c --dask-num-workers {threads} \
            contacts export {input.contacts} \
            salsa_bed {params.prefix}  2>{log}
        set +u; conda deactivate; set -u
    """

rule salsaBed2bam:
    input:
        salsa = rules.to_salsa_bed.output,
        ref = fa,
        chrsize = chromsizes,
        script = config["pip_dir"] + 'salsaBed2bam.sh',
    output:
        "salsa/MboI_{sample}_{ref}_{phase}_filter.salsa.bam",
    shell: """
        bash {input.script} {input.salsa} \
            {input.ref} {input.chrsize}
    """

rule merge_bam:
    input:
        expand(rules.salsaBed2bam.output,
        sample = samples, ref = reference, phase = phases),
    output:
        "merged_salsa.bam",
    shell: """
        set +u; source activate snakemake; set -u
        samtools merge {output} {input}
        set +u; conda deactivate; set -u
        """

hbf_dir = '/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/database/hic_breakfinder'
rule hic_breakfinder:
    input:
        bam = rules.merge_bam.output,
        inter = hbf_dir + 'inter_expect_1Mb.hg38.txt',
        intra = hbf_dir + 'intra_expect_100kb.hg38.txt',
    output:
        "hic_breakfinder/merged_salsa.breaks.txt",
    params:
        "hic_breakfinder/merged_salsa",
    shell: """ 
        set +u; source activate neoloop; set -u
        hic_breakfinder --bam-file {input} \
        --exp-file-inter $inter \
        --exp-file-intra $intra \
        --name {params}
        set +u; conda deactivate; set -u
    """

rule merge_cool:
    input:
        mcool = expand("cooler/MboI_{sample}_{ref}_{phase}_chronly_sort.mcool",
        sample = samples, ref = reference, phase = phases),
    output:
        "merged.cool",
    shell: """
        set +u; source activate neoloop; set -u
        ls {input.mcool} |
        sed 's/mcool/mcool\\:\\:\\/resolutions\\/{wildcards.res}/g' |
        xargs cooler merge {output}
        set +u; conda deactivate; set -u
    """
    
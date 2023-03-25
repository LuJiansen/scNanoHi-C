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
        expand("merged_{ref}_{phase}_filter.salsa.bed",
        ref = reference, phase = phases),

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

rule merge_bed:
    input:
        expand("MboI_{{sample}}_{ref}_{phase}_filter.salsa.bed",sample = samples),
    output:
        "merged_{ref}_{phase}_filter.salsa.bed",
    shell: """
        cat {input} > {output}
    """
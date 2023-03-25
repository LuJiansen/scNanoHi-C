import os
import sys

sample_list = 'sample'
samples=[]
with open(sample_list,'r') as f:
	for line in f:
		samples.append(line.strip('\n'))
print(samples)
phases = config["phase"]
refs =config["reference"]

wildcard_constraints:
    sample="[0-9A-Za-z]+"

rule all:
    input:
        expand("MboI_{sample}_{ref}_{phase}.mono.txt",
        sample = samples,phase = phases, ref = refs),
        expand("MboI_{sample}_{ref}_{phase}.mono.cov",
        sample = samples,phase = phases, ref = refs),
        expand("cells_{ref}_{phase}_chr_cov_summary.csv",
        phase = phases, ref = refs),
        expand("cells_{ref}_{phase}_genomecov_summary.csv",
        phase = phases, ref = refs),
        "merged_monomer.txt.gz",
        "merged_bincov_summary.txt",
    
rule pq2df:
    input:
        pq = "MboI_{sample}_{ref}_{phase}.contacts.clean.parquet",
        script = config["pip_dir"] + 'pq2monomer.R',
    params:
        prefix = config["prefix"],
    output:
        "MboI_{sample}_{ref}_{phase}.mono.txt",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.sample} \
            {input.pq} \
            {output} \
            {params.prefix}
        set +u; conda deactivate; set -u
    """

rule merge_mono:
    input:
        expand(rules.pq2df.output,sample = samples, phase = phases, ref = refs),
    output:
        "merged_monomer.txt.gz",
    params:
        "merged_monomer.txt",
    shell: """
        echo -e "chrom\\tstart\\tend\\tread_name\\tread_length\\thaplotype\\tread_idx\\tcell" > {params}
        cat {input} >> {params}
        pigz {params}
    """
    
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

rule bincov:
    input:
        mono = rules.pq2df.output,
        bin_1M = config["pip_dir"] + refs + '_1M_windows.bed',
        bin_500k = config["pip_dir"] + refs + '_500k_windows.bed',
        bin_100k = config["pip_dir"] + refs + '_100k_windows.bed',
        bin_50k = config["pip_dir"] + refs + '_50k_windows.bed',
    output:
        temp("MboI_{sample}_{ref}_{phase}_bincov.txt"),
    shell: """
        bincov(){{
            cov=$1
            mono=$2
            n=`bedtools intersect -a $cov -b $mono -wa |
            uniq | wc -l`
            echo $n
        }}

        a=`bincov {input.bin_1M} {input.mono}`
        b=`bincov {input.bin_500k} {input.mono}`
        c=`bincov {input.bin_100k} {input.mono}`
        d=`bincov {input.bin_50k} {input.mono}`

        echo -e "{wildcards.sample}\\t$a\\t$b\\t$c\\t$d" > {output}
    """
    
rule collect_bincov:
    input:
        expand(rules.bincov.output,
        sample = samples,phase = phases, ref = refs)
    output:
        "merged_bincov_summary.txt",
    shell: """
        echo -e "sample\\tbin_1M\\tbin_500k\\tbin_100k\\tbin_50k" > {output}
        cat {input} >> {output}
    """
    
rule df2bedcov:
    input:
        mono = rules.pq2df.output,
        chr_size = "../refgenome/{ref}.rg.chromsizes"
    output:
        temp("MboI_{sample}_{ref}_{phase}.mono.bed"),
    shell: """
        cat {input.mono} | 
        awk '{{print $1,$2,$3}}' OFS='\\t' |
        sort -k 1,1 | 
        bedtools genomecov -bg -i stdin -g {input.chr_size} |
        awk '{{print $1,$2,$3,$4,$3-$2,$4*($3-$2)}}' OFS='\\t' > {output}
    """

rule bedcov:
    input:
        rules.df2bedcov.output,
    output:
        "MboI_{sample}_{ref}_{phase}.mono.cov",
    shell: """
        bedtools groupby -i {input} \
            -g 1 -c 5,6 -o sum,sum > {output}
    """

rule collect_bedcov:
    input:
        dep = expand("MboI_{sample}_{{ref}}_{{phase}}.mono.cov",sample = samples),
        chrSize = "../refgenome/{ref}.rg.chromsizes",
        script = config["pip_dir"] + 'collect_coverage.R',
    output:
        chr_cov = "cells_{ref}_{phase}_chr_cov_summary.csv",
        genomecov = "cells_{ref}_{phase}_genomecov_summary.csv",
    params:
        chrregx = "'chr[0-9]+|[XY]$'",
        cov_path = '.',
    shell: """ 
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {params.cov_path} \
            {input.chrSize} \
            {wildcards.ref} \
            {wildcards.phase} \
            {params.chrregx}
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
    
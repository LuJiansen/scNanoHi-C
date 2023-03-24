sample_list = 'sample'
samples=[]
with open(sample_list,'r') as f:
	for line in f:
		samples.append(line.strip('\n'))
print(samples)

pip_dir="/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/Github/Pipeline/scHic/"
phases = "unphased"
reference = "GRCh38"
resolutions = ['100000']

rule all:
    input:
        expand("{sample}_{ref}_{phase}_{res}_adjnTIF.rds",
        sample = samples, ref = reference, phase = phases, res = resolutions),

rule cool2mtx:
    input:
        "cooler/MboI_{sample}_{ref}_{phase}_chronly_sort.mcool",
    output:
        "cooler/MboI_{sample}_{ref}_{phase}_chronly_sort_{res}_all_mtx.gz",
    shell: """
        cool2matrix {input} all {wildcards.res}
        """

rule cool2bins:
    input:
        expand("cooler/MboI_{sample}_{ref}_{phase}_chronly_sort.mcool",
        sample = samples, ref = reference, phase = phases),
    output:
        "cooler/cool_{res}_bins.mtx.gz",
    shell: """
        ls {input} | head -n 1 |
        xargs -I {{}} cool2matrix {{}} bins {wildcards.res} {output}
        """

rule adjnTIF:
    input:
        mtx = rules.cool2mtx.output,
        bins = rules.cool2bins.output,
        scripts = pip_dir + 'cal_adjnTIF.R',
        cnv = "cnv/MboI_{sample}_{ref}_{phase}_chronly_{res}_cnv.bg",
    output:
        "{sample}_{ref}_{phase}_{res}_adjnTIF.rds",
    log:
        "logs/{sample}_{ref}_{phase}_{res}_adjnTIF.log",
    params:
        "{sample}_{ref}_{phase}",
    threads: 10,
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.scripts} \
            {input.bins} \
            {input.mtx} \
            {input.cnv} \
            {wildcards.res} \
            {params} \
            {threads} &> {log}
        set +u; conda deactivate; set -u
    """



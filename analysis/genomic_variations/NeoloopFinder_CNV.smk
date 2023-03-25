sample_list = 'sample'
samples=[]
with open(sample_list,'r') as f:
	for line in f:
		samples.append(line.strip('\n'))
print(samples)
reference = 'GRCh38'
phases = 'unphased'
enzymes = 'MboI'
ploidy = 2

resolutions = ['100000','250000','500000','1000000']

base_dir = '/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/'
pip_dir = base_dir + 'Github/Pipeline/scHic/'
neoloop_resource = base_dir + "database/neoloop_data"
# fragments = '../virtual_digest/{enzyme}_{ref}.vd.fragments.parquet'
# chromsizes = '../refgenome/{ref}.rg.chromsizes'

rule all:
    input:
        expand("cooler/{enzyme}_{sample}_{ref}_{phase}_chronly_sort.mcool",
        sample = samples, ref = reference, phase = phases, enzyme = enzymes),
        expand("cnv/{enzyme}_{sample}_{ref}_{phase}_chronly_{res}_cnv.bg",
        sample = samples, ref = reference, phase = phases,res = resolutions, enzyme = enzymes),
        expand("cnv/{enzyme}_{sample}_{ref}_{phase}_chronly_{res}_cnv_seg.txt",
        sample = samples, ref = reference, phase = phases,res = resolutions, enzyme = enzymes),

rule chronly:
    input:
        pq = "{enzyme}_{sample}_{ref}_{phase}.contacts.clean.parquet",
        script = pip_dir + 'chronly.R',
    output:
        temp("{enzyme}_{sample}_{ref}_{phase}.contacts.chronly.parquet"),
    params:
        regex = '"^chr"',
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.chronly.log",
    threads: 5,
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {input.pq} {params.regex} {output} > {log}
        set +u; conda deactivate; set -u
        """

rule fragments_flt:
    input:
        frag = '../virtual_digest/{enzyme}_{ref}.vd.fragments.parquet',
        script = pip_dir + 'frags_flt.R',
    output:
        "{enzyme}_{ref}_fragements_chronly.parquet",
    params:
        regex = '"^chr"',
    shell: """ 
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {input.frag} {params.regex} {output}
        set +u; conda deactivate; set -u
        """

rule chromsize_flt:
    input:
        '../refgenome/{ref}.rg.chromsizes'
    output:
        flt = "{ref}_chromsizes_chronly",
        sort = "{ref}_chromsizes_chronly_sort",
    params:
        regex = '"^chr"',
    shell: """ 
        cat {input} | grep {params.regex} > {output.flt}
        sort -V {output.flt} > {output.sort}
        """

rule to_cooler:
    input:
        contacts = rules.chronly.output,
        fragments = rules.fragments_flt.output,
        chromsizes = rules.chromsize_flt.output.flt,
    output:
        temp("cooler/{enzyme}_{sample}_{ref}_{phase}_chronly.cool"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.tocooler.log",
    params:
        prefix = "cooler/{enzyme}_{sample}_{ref}_{phase}_chronly",
    threads: 5,
    shell: """
        set +u; source activate HiC; set -u
        pore_c --dask-num-workers {threads} \
            contacts export {input.contacts} \
            cooler {params.prefix} \
            --fragment-table {input.fragments} \
            --chromsizes {input.chromsizes} 2>{log}
        set +u; conda deactivate; set -u
        """

rule cool_sort:
    input:
        cool = rules.to_cooler.output,
        chrsize = rules.chromsize_flt.output.sort,
    output:
        temp("cooler/{enzyme}_{sample}_{ref}_{phase}_chronly_sort.cool"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.cool_sort.log",
    shell: """
        set +u; source activate HiC; set -u
        cooler dump --join {input.cool} |
        cooler load --format bg2 {input.chrsize}:1000 - {output}
        set +u; conda deactivate; set -u
    """

rule to_mcool:
    input:
        cool = rules.cool_sort.output,
    output:
        "cooler/{enzyme}_{sample}_{ref}_{phase}_chronly_sort.mcool",
    params:
        "10000,20000,50000,100000,250000,500000,1000000,2500000",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.tomcooler.log",
    shell: """
        set +u; source activate HiC; set -u
        cooler zoomify -r {params} {input} 
        set +u; conda deactivate; set -u
    """

rule calc_cnv:
    input:
        mcool = rules.to_mcool.output,
        cache = neoloop_resource,
    output:
        "cnv/{enzyme}_{sample}_{ref}_{phase}_chronly_{res}_cnv.bg",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_{res}.calc_cnv.log",
    shell: """
        set +u; source activate neoloop; set -u
        if [ {wildcards.ref} == "GRCh38" ];then
            genome=hg38
        elif [ {wildcards.ref} == "GRCh37" ];then
            genome=hg37
        elif [ {wildcards.ref} == "GRCm37" ];then
            genome=mm9
        elif [ {wildcards.ref} == "GRCm38" ];then
            genome=mm10
        else
            genome={wildcards.ref}
        fi

        calculate-cnv -H {input.mcool}::resolutions/{wildcards.res} \
        --output {output} \
        -g $genome \
        -e {wildcards.enzyme} \
        --cachefolder {input.cache} \
        --logFile {log}
        set +u; conda deactivate; set -u
    """

rule seg_cnv:
    input:
        rules.calc_cnv.output,
    output:
        "cnv/{enzyme}_{sample}_{ref}_{phase}_chronly_{res}_cnv_seg.txt",
    params:
        ploidy,
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_{res}.seg_cnv.log",
    shell: """
        set +u; source activate neoloop; set -u
        segment-cnv --cnv-file {input} \
        --output {output} \
        --binsize {wildcards.res} \
        --logFile {log} \
        --ploidy {params}
        set +u; conda deactivate; set -u
    """
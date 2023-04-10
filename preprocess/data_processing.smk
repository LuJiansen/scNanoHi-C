pip_dir = "/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/Github/Pipeline/scHic/"
hic_env = "hic_env"

SAMPLE = None # DO NOT modified this

def raw_data_files(wildcards):
    ck_output = checkpoints.demultiplex.get(**wildcards).output[0]
    global SAMPLE
    SAMPLE, = glob_wildcards(os.path.join(ck_output, "{sample}.fastq.gz"))
    return expand(os.path.join(ck_output, "{SMP}.fastq.gz"),SMP=SAMPLE)

def trim_data_files(wildcards):
    ck_output = checkpoints.demultiplex.get(**wildcards).output[0]
    SAMPLE, = glob_wildcards(os.path.join(ck_output, "{sample}.fastq.gz"))
    return expand(os.path.join("trim", "{SMP}_trimed.fastq.gz"), SMP=SAMPLE)

rule all:
    input:
        raw_data_files,
        trim_data_files,
        "raw_fastq.stats",
        "trim_fastq.stats",
        "datasize_summary",
        
rule barcode_generate:
    input:
        index = "bc_index",
        script = pip_dir + 'barcode.sh',
    output:
        "barcode.fa",
    shell: """
        sh {input.script} {input.index}
        """

checkpoint demultiplex:
    input:
        barcode = rules.barcode_generate.output,
        script = pip_dir + 'nanoplexer.sh',
        fq = "pass.fastq.gz",
    output:
        directory("raw_data"),
    params:
        outdir = "raw_data",
        threads_used = 16,
        env = hic_env,
    threads: 20,
    shell: """
        set +u; source activate {params.env}; set -u
        nanoplexer \
            -b {input.barcode} \
            -t {params.threads_used} \
            -p {params.outdir} \
            {input.fq}
        pigz {params.outdir}/*.fastq
        set +u; conda deactivate; set -u
    """

rule cutadapt:
    input:
        "raw_data/{sample}.fastq.gz",
    output:
        "trim/{sample}_trimed.fastq.gz",
    log:
        "log/{sample}_cutadapt.log",
    threads: 5,
    params:
        env = hic_env,
    shell: """
        set +u; source activate {params.env}; set -u
        cutadapt \
        -g AATGATACGGCGACCACCGAGATCT \
        -g TCGTCGGCAGCGTC \
        -g AGATGTGTATAAGAGACAG \
        -a AGATCTCGGTGGTCGCCGTATCATT \
        -a GACGCTGCCGACGA \
        -a CTGTCTCTTATACACATCT \
        --times 2 \
        --minimum-length 500 \
        -e 0.2 \
        -O 12 \
        -j 4 \
        -o {output} {input} \
        > {log}
        set +u; conda deactivate; set -u
    """

# this rule reports error in snakemake=5.5.4; work in 5.26.1
rule raw_data_stats:
    input:
        fq = raw_data_files,
    output:
        "raw_fastq.stats",
    threads: 10,
    params:
        env = hic_env,
    shell: """
        set +u; source activate {params.env}; set -u
        seqkit stats -a -T -j {threads} \
            {input.fq} > {output}
        set +u; conda deactivate; set -u
    """

rule trim_data_stats:
    input:
        fq = trim_data_files,
    output:
        "trim_fastq.stats",
    threads: 10,
    params:
        env = hic_env,
    shell: """
        set +u; source activate {params.env}; set -u
        seqkit stats -a -T -j {threads} \
            {input.fq} > {output}
        set +u; conda deactivate; set -u
    """

rule check_datasize:
    input:
        raw = rules.raw_data_stats.output,
        trim = rules.trim_data_stats.output,
    output:
        "datasize_summary",
    shell: """
        raw_used=`cat raw_fastq.stats | sed 1d |
        grep -v un | awk '{{sum+=$5}};END{{print sum}}'`
        raw_total=`cat raw_fastq.stats | sed 1d |
        awk '{{sum+=$5}};END{{print sum}}'`
        raw_un=`cat raw_fastq.stats | sed 1d | grep un |
        awk '{{sum+=$5}};END{{print sum}}'`

        trim_total=`cat trim_fastq.stats | sed 1d |
        awk '{{sum+=$5}};END{{print sum}}'`
        trim_used=`cat trim_fastq.stats | sed 1d |
        grep -v un | awk '{{sum+=$5}};END{{print sum}}'`
        trim_un=`cat trim_fastq.stats | sed 1d | grep un |
        awk '{{sum+=$5}};END{{print sum}}'`

        echo -e "raw_total\\t$raw_total" > {output}
        echo -e "raw_classified\\t$raw_used" >> {output}
        echo -e "raw_unclassified\\t$raw_un" >> {output}
        echo -e "trim_total\\t$trim_total" >> {output}
        echo -e "trim_classified\\t$trim_used" >> {output}
        echo -e "trim_unclassified\\t$trim_un" >> {output}
    """
import os
import sys

sample_list = 'sample'
samples=[]
with open(sample_list,'r') as f:
	for line in f:
		samples.append(line.strip('\n'))
print(samples)
enzyme = config["enzyme"]
phases = config["phase"]
reference = config["reference"]
reps = config['reps']
resolutions = config['resolution']
run_dipc = config['run_dipc']
run_model = config['run_model']

wildcard_constraints:
    sample="[0-9A-Za-z]+"

print(enzyme)
print(reference)
print(phases)
print(reps)
print(resolutions)
print(run_dipc)
# need wildcard for enzymes

fragments = '../virtual_digest/' + enzyme + '_' + reference + '.vd.fragments.parquet'
chromsizes = '../refgenome/' + reference + '.rg.chromsizes'
print(fragments)
print(chromsizes)

# filtration output files
filter_pq = expand("{enzyme}_{sample}_{ref}_{phase}.contacts.filter.parquet",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),
clean_pq = expand("{enzyme}_{sample}_{ref}_{phase}.contacts.clean.parquet",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),

filter_cool = expand("cooler/{enzyme}_{sample}_{ref}_{phase}_filter.cool",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),
high_cool = expand("cooler/{enzyme}_{sample}_{ref}_{phase}_high.cool",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),
high_unflt_cool = expand("cooler/{enzyme}_{sample}_{ref}_{phase}_high_unflt.cool",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),

# summary
stats = expand("merged_{enzyme}_{ref}_{phase}_stats_summary.csv",phase = phases,ref = reference, enzyme = enzyme),

# dip-c 2D outputs
filter_pairs = expand("pairs/{enzyme}_{sample}_{ref}_{phase}_filter_sort.pairs.gz",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),
cpg_2d = expand("dip_c/{enzyme}_{sample}_{ref}_{phase}_filter_sort.cpg_b1m.color2",
    sample = samples,phase = phases,ref = reference, enzyme = enzyme),
cpg_2d_sum = expand("cpg_2d_{enzyme}_{ref}_{phase}_b1m_summary.csv",
    phase = phases,ref = reference, enzyme = enzyme),

# modelling output files
model_align = expand("dip_c/{enzyme}_{sample}_{ref}_{phase}.{resolution}.align.color",
    sample = samples,phase = phases,resolution = resolutions,ref = reference, enzyme = enzyme),
chr_cif = expand("dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.n.cif",
    sample = samples,phase = phases,rep = reps,resolution = resolutions,ref = reference, enzyme = enzyme),

RMSD_summary = "RMSD_summary.txt",
cpg_3d_sum = expand("cpg_3d_{enzyme}_{ref}_{phase}_{resolution}_summary.csv",
    phase = phases,resolution = resolutions,ref = reference, enzyme = enzyme),
radi_pos_summary = expand("radical_pos_{enzyme}_{ref}_{phase}_{resolution}_summary.csv",
            phase = phases,resolution = resolutions,ref = reference, enzyme = enzyme),
intra_perc_summary = expand("intra_perc_{enzyme}_{ref}_{phase}_{resolution}_summary.csv",
            phase = phases,resolution = resolutions,ref = reference, enzyme = enzyme),

filter_output = [filter_pq,clean_pq,filter_cool,high_cool,high_unflt_cool,stats]
dipc_2d_output = [filter_pairs,stats,cpg_2d,cpg_2d_sum]
model_output = [model_align,chr_cif,RMSD_summary,cpg_3d_sum,radi_pos_summary,intra_perc_summary]

if run_dipc:
    rule_all_output = filter_output + dipc_2d_output
    print("run dip-c")
    if phases != "unphased" and run_model:
        rule_all_output.extend(model_output)
        print("perform modelling")
else:
    rule_all_output = filter_output
    print("filtration only")

#print(rule_all_output)

rule all:
    input:
        data = rule_all_output,

rule filter:
    input:
        pq = "{enzyme}_{sample}_{ref}_{phase}.contacts.parquet",
        script = config["pip_dir"] + 'contact_filtration.R',
    output:
        "{enzyme}_{sample}_{ref}_{phase}.contacts.filter.parquet",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.contacts.log",
    threads: 5,
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.sample} {input.pq} {output} > {log}
        set +u; conda deactivate; set -u
        """

rule apply_filter:
    input:
        pq = rules.filter.output,
        script = config["pip_dir"] + 'apply_filter.R',
    output:
        "{enzyme}_{sample}_{ref}_{phase}.contacts.clean.parquet",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.applyFlt.log",
    threads: 5,
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {input.pq} {output} > {log}
        set +u; conda deactivate; set -u
        """

rule to_cooler:
    input:
        contacts = rules.apply_filter.output,
        fragments = fragments,
        chromsizes = chromsizes,
    output:
        "cooler/{enzyme}_{sample}_{ref}_{phase}_filter.cool",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.cooler.log",
    params:
        prefix = "cooler/{enzyme}_{sample}_{ref}_{phase}_filter",
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

rule high_order:
    input:
        unflt = "{enzyme}_{sample}_{ref}_{phase}.contacts.parquet",
        flt = rules.apply_filter.output,
        script = config["pip_dir"] + 'high_order_contacts.R',
    output:
        flt = temp("{enzyme}_{sample}_{ref}_{phase}.contacts.high.parquet"),
        unflt = temp("{enzyme}_{sample}_{ref}_{phase}.contacts.unflt.high.parquet"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.high_order.log",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {input.flt} {output.flt} > {log}
        Rscript {input.script} \
            {input.unflt} {output.unflt} >> {log}
        set +u; conda deactivate; set -u
        """

rule high_order_cool:
    input:
        flt = rules.high_order.output.flt,
        unflt = rules.high_order.output.unflt,
        fragments = fragments,
        chromsizes = chromsizes,
    output:
        flt = "cooler/{enzyme}_{sample}_{ref}_{phase}_high.cool",
        unflt = "cooler/{enzyme}_{sample}_{ref}_{phase}_high_unflt.cool",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.high_cooler.log",
    params:
        prefix_flt = "cooler/{enzyme}_{sample}_{ref}_{phase}_high",
        prefix_unflt = "cooler/{enzyme}_{sample}_{ref}_{phase}_high_unflt",
    threads: 5,
    shell: """
        set +u; source activate HiC; set -u
        pore_c --dask-num-workers {threads} \
            contacts export {input.flt} \
            cooler {params.prefix_flt} \
            --fragment-table {input.fragments} \
            --chromsizes {input.chromsizes} 2>{log}
        pore_c --dask-num-workers {threads} \
            contacts export {input.unflt} \
            cooler {params.prefix_unflt} \
            --fragment-table {input.fragments} \
            --chromsizes {input.chromsizes} 2>>{log}
        set +u; conda deactivate; set -u
        """

rule filter_stat:
    input:
        pq = rules.filter.output,
        script = config["pip_dir"] + 'filter_stat.R',
    output:
        stat = temp("{sample}_{enzyme}_{ref}_{phase}_contact_filtration_stats.txt"),
        high_order = temp("{sample}_{enzyme}_{ref}_{phase}_order_stats.txt"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.filter.log",
    threads: 2,
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.sample} \
            {wildcards.enzyme} \
            {wildcards.ref} \
            {wildcards.phase} \
            {input.pq} > {log}
        set +u; conda deactivate; set -u
        """

rule porec_stat:
    input:
        mapping_dir = directory("../mapping"),
        align_table_dir = directory("../align_table"),
        cont = "{enzyme}_{sample}_{ref}_{phase}.concatemer_summary.csv",
        script = config["pip_dir"] + 'porec_stat.R',
    output:
        temp("{sample}_{enzyme}_{ref}_{phase}_pore_c_stats.txt"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.porec_stat.log",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.sample} \
            {wildcards.enzyme} \
            {wildcards.ref} \
            {wildcards.phase} \
            {input.cont} > {log}
        set +u; conda deactivate; set -u
    """

rule merge_stat:
    input:
        stat1 = expand("{sample}_{{enzyme}}_{{ref}}_{{phase}}_contact_filtration_stats.txt",sample = samples),
        stat2 = expand("{sample}_{{enzyme}}_{{ref}}_{{phase}}_pore_c_stats.txt",sample = samples),
        high_order = expand("{sample}_{{enzyme}}_{{ref}}_{{phase}}_order_stats.txt",sample = samples),
        # stat1 = expand(rules.filter_stat.output,sample = samples,allow_missing=True),
        # stat2 = expand(rules.porec_stat.output,sample = samples,allow_missing=True),
        script = config["pip_dir"] +'merge_stat.R',
    output:
        all_stat = "merged_{enzyme}_{ref}_{phase}_stats_all.csv",
        selected = "merged_{enzyme}_{ref}_{phase}_stats_selected.csv",
        sum_stat = "merged_{enzyme}_{ref}_{phase}_stats_summary.csv",
        high_order = "merged_{enzyme}_{ref}_{phase}_high_order_summary.csv",
    params:
        phase = config['phase'],
    log:
        "logs/merged_{enzyme}_{ref}_{phase}_stat.log",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
                {wildcards.enzyme} \
                {wildcards.ref} \
                {wildcards.phase} \
                2> {log}
        set +u; conda deactivate; set -u
    """

rule to_phased_pairs:
    input:
        contacts = rules.apply_filter.output,
        script = config["pip_dir"] + 'contacts_to_phased_pairs.R',
        chromsizes = chromsizes,
    output:
        temp("pairs/{enzyme}_{sample}_{ref}_{phase}_filter.pairs"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.to_pairs.log",
    params:
        tmp = "pairs/{enzyme}_{sample}_{ref}_{phase}_filter.tmp",
        header = "pairs/{enzyme}_{sample}_{ref}_{phase}_filter.pairs.header"
    threads: 5,
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {input.contacts} \
            {params.tmp} > {log}
        echo "## pairs format v1.0" > {params.header}
        echo "#shape: upper triangle" >> {params.header}
        sed 's/\t/ /g' {input.chromsizes} |
            sed 's/^/\#chromosome\: /g' >> {params.header}
        echo "#columns: readID chr1 pos1 chr2 pos2 strand1 strand2 phase0 phase1" >> {params.header}
        
        cat {params.header} {params.tmp} > {output}
        rm {params.header} {params.tmp}
        set +u; conda deactivate; set -u
    """

rule sort_pairs:
    input:
        rules.to_phased_pairs.output,
    output:
        temp("pairs/{enzyme}_{sample}_{ref}_{phase}_sort.pairs.gz"),
    threads: 2,
    shell: """
        set +u; source activate HiC; set -u
        pairtools sort {input} --nproc {threads} | pigz > {output}
        set +u; conda deactivate; set -u
    """

# need option of gender and par_bed
# will cause error in mixed genome which seqnames start with not chr
rule filter_pairs:
    input:
        pair = rules.sort_pairs.output,
        par_bed = config['par_bed'],
    output:
        "pairs/{enzyme}_{sample}_{ref}_{phase}_filter_sort.pairs.gz",
    params:
        gender = config['gender'],
    shell: """
        set +u; source activate LJS_py2; set -u
        # for female
        if [ {params.gender} == "female" ];then  
            hickit.js chronly -y {input.pair} | pigz > {output}
        elif [ {params.gender} == "male" ];then
            # for male
            hickit.js chronly {input.pair} | hickit.js bedflt {input.par_bed} - | pigz > {output}
        fi
        set +u; conda deactivate; set -u
    """

rule hickit_impute:
    input:
        pairs = rules.filter_pairs.output,
        hickit = config['hickit'],
    output:
        "hickit/{enzyme}_{sample}_{ref}_{phase}.imput.pairs.gz",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}.hickit_impute.log",
    params:
        "hickit/{enzyme}_{sample}_{ref}_{phase}.imput.pairs",
    shell: """
        {input.hickit} -i {input.pairs} -u -o {params} 2>{log}
        pigz {params}
    """

rule hickit_model:
    input:
        pairs = rules.hickit_impute.output,
        hickit = config['hickit'],
    output:
        tdg1 = "hickit/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.1m.3dg",
        tdg2 = "hickit/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.100k.3dg",
        tdg3 = "hickit/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.50k.3dg",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_rep{rep}.hickit_model.log",
    shell: """
        {input.hickit} -M -s{wildcards.rep} -i {input.pairs} \
            -Sr1m -c1 -r10m -c2 \
            -b4m -b1m -O {output.tdg1} \
            -D5 -b100k -O {output.tdg2} \
            -D5 -b50k -O {output.tdg3} 2>>{log}
    """

rule pair2dipc:
    input:
        pairs = rules.filter_pairs.output,
        pair2con = config['dipc_dir'] + 'scripts/hickit_pairs_to_con.sh',
    output:
        pairs = "pairs/{enzyme}_{sample}_{ref}_{phase}_filter_sort.con.gz",
    shell: """
        {input.pair2con} {input.pairs}
    """

rule impute2dipc:
    input:
        impute = rules.hickit_impute.output,
        impute2con = config['dipc_dir'] + 'scripts/hickit_impute_pairs_to_con.sh',
    output:
        impute = "hickit/{enzyme}_{sample}_{ref}_{phase}.imput.con.gz",
    shell: """
        {input.impute2con} {input.impute}
    """

rule clean_cons:
    input:
        tdg = "hickit/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.3dg",
        rescale = config['dipc_dir'] + 'scripts/hickit_3dg_to_3dg_rescale_unit.sh',
        impute = rules.impute2dipc.output.impute,
    output:
        tmp = temp("hickit/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.dip-c.3dg"),
        tdg = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.clean.3dg",
    shell: """
        set +u; source activate LJS_py2; set -u
        {input.rescale} {input.tdg}
        dip-c clean3 -c {input.impute} {output.tmp} > {output.tdg}
        set +u; conda deactivate; set -u
    """

rule dipc_align:
    input:
        expand( "dip_c/{{enzyme}}_{{sample}}_{{ref}}_{{phase}}.rep{rep}.{{resolution}}.clean.3dg",rep = reps),
        # expand(rules.clean_cons.output.tdg,rep = reps,allow_missing=True),
    output:
        "dip_c/{enzyme}_{sample}_{ref}_{phase}.{resolution}.align.color",
    params:
        "dip_c/{enzyme}_{sample}_{ref}_{phase}.{resolution}.",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_{resolution}_dipc_align.log",
    shell: """
        set +u; source activate LJS_py2; set -u
        dip-c align -o {params} \
            {input} \
            2> {log} \
            > {output}
        set +u; conda deactivate; set -u
    """

rule get_RMSD:
    input:
        dep =  rules.dipc_align.output,
        log = rules.dipc_align.log,
    output:
        temp("dip_c/{enzyme}_{sample}_{ref}_{phase}.{resolution}.RMSD.tmp"),
    shell: """
        paste <(echo "{wildcards.sample}_{wildcards.phase}_{wildcards.resolution}") \
        <(cat {input.log} | grep median | sed 's/.*: //g' |
        tr "\\n" "\\t" | sed 's/\\t$//g') \
        <(cat {input.log} | grep "RMS " | sed 's/.*: //g') > {output}
    """

rule RMSD_summary:
    input:
        expand(rules.get_RMSD.output,sample = samples,
        phase = phases,resolution = resolutions,ref=reference, enzyme = enzyme),
    output:
        "RMSD_summary.txt",
    shell: """
        echo -e "sample\\tround1\\tround2\\tround3\\tmedian\\tRMS" > {output}
        cat {input} >> {output}
    """

rule dipc_2d_color:
    input:
        contact = rules.pair2dipc.output.pairs,
        cpg1m = config["dipc_dir"] + 'color/{ref}.cpg.1m.txt',
    output:
        "dip_c/{enzyme}_{sample}_{ref}_{phase}_filter_sort.cpg_b1m.color2",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_dipc_2d_color.log",
    shell: """
        set +u; source activate LJS_py2; set -u
        dip-c color2 -b1000000 -H -c {input.cpg1m} \
            -s {input.contact} > {output} 2> {log}
        set +u; conda deactivate; set -u
    """

rule cpg_2d_summary:
    input:
        cpg = expand( "dip_c/{{enzyme}}_{sample}_{{ref}}_{{phase}}_filter_sort.cpg_b1m.color2",sample = samples),
        # cpg = expand(rules.dipc_2d_color.output,sample = samples,allow_missing = True),
        script = config['pip_dir'] + 'color_score.R',
    output:
        "cpg_2d_{enzyme}_{ref}_{phase}_b1m_summary.csv",
    log:
        "logs/cpg_2d_{enzyme}_{ref}_{phase}_summary.log",
    params:
        indir = "dip_c",
        pattern1 = "_{ref}_{phase}",
        pattern2 = "_filter_sort.cpg_b1m.color2",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.enzyme} \
            {params.indir} \
            {params.pattern1} \
            {params.pattern2} \
            {output} > {log}
        set +u; conda deactivate; set -u
    """

rule dipc_3d_color:
    input:
        tdg = rules.clean_cons.output.tdg,
        cpg = config["dipc_dir"] + 'color/{ref}.cpg.{resolution}.txt',
        chrfile = config["dipc_dir"] + 'color/{ref}.chr.txt',
    output:
        cpg_color = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.cpg_s3.color",
        cpg_cif = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.cpg_s3.cif",
        n_cif = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.n.cif",
        exp = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.exp_n.3dg",
        exp_cif = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.exp_n.cif",
        label = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.label_n.3dg",
        label_cif = "dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.label_n.cif",
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_rep{rep}_{resolution}_dipc_3d_color.log",
    shell: """
        set +u; source activate LJS_py2; set -u
        dip-c color -c {input.cpg} -s3 {input.tdg} > {output.cpg_color} 2> {log}
        dip-c vis -c {output.cpg_color} {input.tdg} > {output.cpg_cif} 2>> {log}
        
        dip-c color -n {input.chrfile} {input.tdg} 2>> {log} |
        dip-c vis -c /dev/stdin {input.tdg} > {output.n_cif} 2>> {log}

        dip-c exp {input.tdg} > {output.exp} 2>> {log}
        dip-c color -n {input.chrfile} {output.exp} 2>> {log} |
        dip-c vis -c /dev/stdin {output.exp} > {output.exp_cif} 2>> {log}

        dip-c exp -c {input.tdg} > {output.label} 2>> {log}
        dip-c vis {output.label} | sed 's/(mat)/♀/g; s/(pat)/♂/g' > {output.label_cif} 2>> {log}
        set +u; conda deactivate; set -u
    """

rule cpg_3d_summary:
    input:
        cpg = expand("dip_c/{{enzyme}}_{sample}_{{ref}}_{{phase}}.rep{rep}.{{resolution}}.cpg_s3.color",sample = samples,rep = reps),
        # cpg = expand(rules.dipc_3d_color.output.cpg_color,sample = samples,rep = reps,allow_missing = True),
        script = config['pip_dir'] + 'color_score.R',
    output:
        "cpg_3d_{enzyme}_{ref}_{phase}_{resolution}_summary.csv",
    log:
        "logs/cpg_3d_{enzyme}_{ref}_{phase}_{resolution}_summary.log",
    params:
        indir = "dip_c",
        pattern1 = "_{ref}_{phase}.",
        pattern2 = ".{resolution}.cpg_s3.color",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.enzyme} \
            {params.indir} \
            {params.pattern1} \
            {params.pattern2} \
            {output} > {log}
        set +u; conda deactivate; set -u
    """

rule dipc_radical_pos:
    input:
        tdg = rules.clean_cons.output.tdg,
    output:
        temp("dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.radical_pos.color"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_rep{rep}_{resolution}_radical_pos.log",
    shell: """
        set +u; source activate LJS_py2; set -u
        dip-c color -C {input.tdg} > {output} 2> {log}
        set +u; conda deactivate; set -u
    """

rule radical_3d_summary:
    input:
        radical = expand("dip_c/{{enzyme}}_{sample}_{ref}_{{phase}}.rep{rep}.{{resolution}}.radical_pos.color",
            sample = samples,rep = reps,ref=reference),
        # radical = expand(rules.dipc_radical_pos.output,
        #     sample = samples,rep = reps,ref=reference,allow_missing = True),
        script = config['pip_dir'] + 'color_score.R',
    output:
        "radical_pos_{enzyme}_{ref}_{phase}_{resolution}_summary.csv",
    log:
        "logs/radical_pos_{enzyme}_{ref}_{phase}_{resolution}_summary.log",
    params:
        indir = "dip_c",
        pattern1 = "_{ref}_{phase}.",
        pattern2 = ".{resolution}.radical_pos.color",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.enzyme} \
            {params.indir} \
            {params.pattern1} \
            {params.pattern2} \
            {output} > {log}
        set +u; conda deactivate; set -u
    """

rule dipc_intra_perc:
    input:
        tdg = rules.clean_cons.output.tdg,
    output:
        temp("dip_c/{enzyme}_{sample}_{ref}_{phase}.rep{rep}.{resolution}.intra.color"),
    log:
        "logs/{enzyme}_{sample}_{ref}_{phase}_rep{rep}_{resolution}_intra.log",
    shell: """
        set +u; source activate LJS_py2; set -u
        dip-c color -i3 {input.tdg} > {output} 2> {log}
        set +u; conda deactivate; set -u
    """

rule intra_perc_summary:
    input:
        radical = expand("dip_c/{{enzyme}}_{sample}_{ref}_{{phase}}.rep{rep}.{{resolution}}.intra.color",
            sample = samples,rep = reps,ref=reference),
        # radical = expand(rules.dipc_intra_perc.output,
        #     sample = samples,rep = reps,ref=reference,allow_missing = True),
        script = config['pip_dir'] + 'color_score.R',
    output:
        "intra_perc_{enzyme}_{ref}_{phase}_{resolution}_summary.csv",
    log:
        "logs/intra_perc_{enzyme}_{ref}_{phase}_{resolution}_summary.log",
    params:
        indir = "dip_c",
        pattern1 = "_{ref}_{phase}.",
        pattern2 = ".{resolution}.intra.color",
    shell: """
        set +u; source activate r_env; set -u
        Rscript {input.script} \
            {wildcards.enzyme} \
            {params.indir} \
            {params.pattern1} \
            {params.pattern2} \
            {output} > {log}
        set +u; conda deactivate; set -u
    """
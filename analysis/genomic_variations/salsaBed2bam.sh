#!/bin/bash
salsa_bed=$1
fa=$2
chrsize=$3

tmpdir=`mktemp -d --tmpdir=.`

# salsa bed to paired bam (without SEQ and QUAL)
cat ${salsa_bed} |
    xargs -n 12 | 
    awk '{if($1==$7){r="=";len=$8-$2}
    else{r=$7;len=0}};
    {if($6=="-"){
        if($12=="-"){f=115}else{f=83}}
    else{
        if($12=="-"){f=99}else{f=67}}
    };{print $4,f,$1,$2+1,$5,$3-$2"M",r,$8+1,len}' OFS=',' > \
    ${tmpdir}/R1.tmp

cat ${salsa_bed} |
    xargs -n 12 | 
    awk '{if($1==$7){r="=";len=$2-$8}
    else{r=$1;len=0}};
    {if($6=="-"){
        if($12=="-"){f=179}else{f=163}}
    else{
        if($12=="-"){f=147}else{f=131}}
    };{print $10,f,$7,$8+1,$11,$9-$8"M",r,$2+1,len}' OFS=',' > \
    ${tmpdir}/R2.tmp

paste ${tmpdir}/R1.tmp ${tmpdir}/R2.tmp |tr "\t" "\n" | tr "," "\t" > ${tmpdir}/salsa.bam.tmp

# bed to SEQ and QUAL dummy
bedtools getfasta -tab -name \
    -fi ${fa} \
    -bed ${salsa_bed}|
awk '{print $1,$2,$2}' |
awk 'gsub('/[A-Za-z]/',"F",$3)' OFS='\t' > \
${tmpdir}/test.salsa.fa

# merged 
paste <(sed "s/\/[12]//g" ${tmpdir}/salsa.bam.tmp) \
    <(awk '{print $2,$3}' OFS='\t' ${tmpdir}/test.salsa.fa) > \
    ${tmpdir}/test.salsa.dummy.sam

# dummy sam header
echo -e "@HD\tVN:1.5\tSO:unknown" > ${tmpdir}/dummy.header
awk '{print "@SQ","SN:"$1,"LN:"$2}' OFS='\t' \
	${chrsize} >> ${tmpdir}/dummy.header

# sam to bam
cat ${tmpdir}/dummy.header ${tmpdir}/test.salsa.dummy.sam | 
samtools view -hb > ${tmpdir}/test.salsa.dummy.bam

mv ${tmpdir}/test.salsa.dummy.bam ${salsa_bed%%bed}bam
rm -r ${tmpdir}
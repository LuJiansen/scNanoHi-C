#!/bin/bash
input=$1
# input contents:
# column1: PCR_barcode
# column2: start TN5_barcode
# column3: end TN5_barcode
# column4: cell name prefix (optional)

pipdir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/lujiansen/Github/Pipeline/scHic
PCR=$pipdir/pcr_barcode.txt
TN5=$pipdir/tn5_barcode.txt
seq1=AATGATACGGCGACCACCGAGATCT
seq2=TCGTCGGCAGCGTC
seq3=AGATGTGTATAAGAGACAG

if [ -s barcode.fa ];then
	echo "Existed barcode.fa was removed!"
	rm barcode.fa
fi

while read i
do
    a=($i)
    if [ ${a[0]} -ge 0 ];then
        PCR_BC=${a[0]}
        BC_start=${a[1]}
        BC_end=${a[2]}
        PCR_Name=${a[3]}
        [ $BC_start ] || BC_start=1
        [ $BC_end ] || BC_end=24
        [ $PCR_Name ] || PCR_Name=P${PCR_BC}
        PCR_Seq=`cat $PCR | awk -v i=${PCR_BC} '$1==i {print $2}'`
        for TN5_BC in `seq ${BC_start} ${BC_end}`
        do
            TN5_Seq=`cat $TN5 | awk -v i=$TN5_BC '$1==i {print $2}'`
            Name=${PCR_Name}B${TN5_BC}
            BARCODE=${seq1}${PCR_Seq}${seq2}${TN5_Seq}${seq3}
            echo "Name: $Name"
            echo "PCR_Seq: $PCR_Seq"
            echo "TN5_Seq: $TN5_Seq"
            echo "=============================================================="
            echo ">${Name}" >> barcode.fa
            echo "$BARCODE" >> barcode.fa
        done
    else
        echo "Barcode is no integer!"
        exit 1
    fi
done < $input

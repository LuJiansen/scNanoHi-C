#!/bin/bash
input=$1
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

for i in `cat $input | cut -f1`
do
	if [ $i -ge 0 ];then
		PCR_BC=$i
		PCR_Name=`cat $input | awk -v i=$i '$1==i {print $2}'`
		[ $PCR_Name ] || PCR_Name=P$i
		PCR_Seq=`cat $PCR | awk -v i=$i '$1==i {print $2}'`
		for TN5_BC in `seq 1 24`
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
done

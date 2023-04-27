#!/usr/bin/sh

dir="/GEO/GSE148822/fastq_files"
sample=`cat $1`

for sam in $sample
	do
		# Step 2: Identify correct cell barcodes
		umi_tools whitelist --stdin $dir/$sam/$sam\_R1.fq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --set-cell-number=100 --log2stderr > $sam\_whitelist.txt

		# Step 3: Extract barcdoes and UMIs and add to read names
		umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin $dir/$sam/$sam\_R1.fq.gz --stdout $sam\_R1_extracted.fastq.gz --read2-in $dir/$sam/$sam\_R2.fq.gz --read2-out=$sam\_R2_extracted.fastq.gz --whitelist=$sam\_whitelist.txt

		echo -e "$sam finished ~~~~~"
done

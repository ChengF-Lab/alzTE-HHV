#!/usr/bin/sh

dir="/GEO/GSE148822"
sample=`cat $1`   ##### This is sample list, each line is file basename

for sam in $sample
	do
		scTE -i $dir/$sam\.bam.1 -o $sam --min_genes 200 --min_counts 400 -p 1 -x /scTE/hg38/hg38.exclusive.idx --keeptmp True
		gzip $sam.csv
		echo -e "$sam finished ~~~~~"
done

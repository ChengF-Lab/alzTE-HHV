#!/bin/bash

samples=`cat $1`  ##### All of your samples, each line is a sample name
input_dir=/Clean_data   ##### The fastq file dir


for sam in $samples
	do
		echo "VIRTUS.SE start at date: "`date`
		mkdir $sam
		cd $sam
		cwltool --tmp-outdir-prefix=${PWD}/tmp_cwl/ --tmpdir-prefix=${PWD}/tmp_cwl/ --rm-tmpdir /VIRTUS2/bin/VIRTUS.SE.cwl --fastq $input_dir/$sam\.fq.gz --genomeDir_human /VIRTUS2/indexs/STAR_index_human --genomeDir_virus /VIRTUS2/indexs/STAR_index_virus --outFileNamePrefix_human human.$sam --filename_output VIRTUS.$sam\.output.tsv
		rm *.bam *.fq
		rm -r tmp_cwl
		cd ../
		echo "$sam finished at date: "`date`
done


echo Successfully completed all samples!



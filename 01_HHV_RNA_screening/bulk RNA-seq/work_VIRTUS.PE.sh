#!/bin/bash

samples=`cat $1`   ###### This file contains all the samples and each line is a sample name
input_dir=/Clean_data  ##### This dir contains all the fastq files


for sam in $samples
	do
		echo "VIRTUS.SE start at date: "`date`
		mkdir $sam
		cd $sam
		cwltool --tmp-outdir-prefix=${PWD}/tmp_cwl/ --tmpdir-prefix=${PWD}/tmp_cwl/ --rm-tmpdir /VIRTUS2/bin/VIRTUS.PE.cwl --fastq1 $input_dir/$sam\.R1.fq.gz --fastq2 $input_dir/$sam\.R2.fq.gz --genomeDir_human /VIRTUS2/indexs/STAR_index_human --genomeDir_virus /VIRTUS2/indexs/STAR_index_virus --outFileNamePrefix_human human.$sam --filename_output VIRTUS.$sam\.output.tsv --nthreads 5
		rm *.bam *.fq
		rm -r tmp_cwl
		cd ../
		echo "$sam finished at date: "`date`
done


echo Successfully completed all samples!



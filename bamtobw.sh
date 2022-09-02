#!/bin/bash
#This script will convert  sorted and indexed BAM files in the working directory to BIGWIG.

# get genome sizes with refgenie
chromSize=$(refgenie seek hg19/fasta.chrom_sizes)

#generate bed files:
for i in *.bam
do 
	# get the BAM file name and clean file name without any extensions
	fileName="${i%.*}"
	BED_name="${fileName}.bed"
	COV_name="${fileName}.cov"
	BW_name="${fileName}.bw"
	
	# get BED coverage files -> cov -> bw
	bedtools bamtobed -i $i > $BED_name
	bedtools genomecov -i $BED_name -g $chromSize -bg > $COV_name
	bedGraphToBigWig $COV_name $chromSize $BW_name
done
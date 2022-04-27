#!/bin/bash

# The script takes BAM files and creates bigWig files and BED coverage files

mkdir output_bw_notNormalized
mkdir output_BED

# get chromosome sizes with refgenie
chromSize=$(refgenie seek hg19/fasta.chrom_sizes)

# go through BAM files - create BED coverage files
# the create unnormalized bw files
for bamFile in *.bam
do 
	# extract file name - everything before the first dot
	fileName="${bamFile%%.*}"
	
	# create file names for outputs
	bedFile="${fileName}.bed"
	covFile="${fileName}.cov"
	bwFile="${fileName}.bw"
	
	# make BED files
	# -> cov files
	#  -> bw files
	bedtools bamtobed -i $bamFile > $bedFile
	bedtools genomecov -i $bedFile -g $chromSize -bg > $covFile; 
	bedGraphToBigWig $covFile $chromSize $bwFile
	
	# remove unnecessary cov file
	rm $covFile
	
	# move files to output directory
	mv $bedFile output_BED
	mv $bwFile output_bw_notNormalized

done







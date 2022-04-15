#!/bin/bash

# the script maps all fastq files within a folder to dm6, removes unmapped reads, and removes blacklisted sites
# outputs are sorted and indexed indexed BAM files


# get bowtie2 index
dm6index=$(refgenie seek dm6/bowtie2_index)
echo -e "Give a full name (including pathway) of a file containing dm6 blacklisted sites."
read blacklistSites


#unzip files:
for zippedFastq in *.fastq.gz 
do 
	# get the fastq file name and clean file name without any extensions
	fastqFile="${zippedFastq%.*}"
	fileName="${zippedFastq%%.*}"
	
	# first unzip the file
	gunzip $zippedFastq
	
	# Map to dm6 with bowtie2
	echo "Mapping following file:"
	echo $fastqFile
	
	unsortedBAM="${fileName}_dm6_unsorted.bam"
	bowtie2 -x $dm6index -p 6 --time -U $fastqFile | samtools view -S -b - > $unsortedBAM
	
	#sort and index the bam file
	sorted="${fileName}_dm6_sorted"
	sortedBAM="${sorted}.bam"
	sortedBAI="${sortedBAM}.bai"
	samtools sort $unsortedBAM $sorted
	samtools index $sortedBAM
	
	# remove the unsorted BAM file
	rm $unsortedBAM

	#filter unmapped reads
	filteredBAM="${fileName}_dm6_filtered.bam"
	echo -e The numbers of mapped and unmapped reads to dm6, respectively:
	samtools view -c -F 4 $sortedBAM
	samtools view -c -f 4 $sortedBAM
	samtools view -h -F 4 -b $sortedBAM > $filteredBAM
	
	# remove unfiltered files
	rm $sortedBAM
	rm $sortedBAI

	#sort and index
	sorted="${fileName}_filtered_sorted"
	sortedBAM="${sorted}.bam"
	sortedBAI="${fileName}_filtered_sorted.bam.bai"
	samtools sort $filteredBAM $sorted
	samtools index $sortedBAM
	
	# remove unsorted file
	rm $filteredBAM

	#filter blacklisted sites & make index for resulting file
	finalBAM="${fileName}_dm6.bam"
	bedtools intersect -abam $sortedBAM -b $blacklistSites -v > $finalBAM
	samtools index $finalBAM
	
	# remove files without blacklist removed
	rm $sortedBAM
	rm $sortedBAI

	echo -e Final Drosophila read count for this dataset after removal of blacklisted sequences:
	samtools view -c -F 4 $finalBAM
	
	# zip the fastq file back
	gzip $fastqFile
done







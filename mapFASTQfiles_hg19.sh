#!/bin/bash

# the script maps all fastq files within a folder to hg19, removes unmapped reads, and removes blacklisted sites
# outputs are sorted and indexed indexed BAM files


# get bowtie2 index
hg19index=$(refgenie seek hg19/bowtie2_index)
echo -e "Give a full name (including pathway) of a file containing hg19 blacklisted sites."
read blacklistSites


#unzip files:
for zippedFastq in *.fastq.gz 
do 
	# get the fastq file name and clean file name without any extensions
	fastqFile="${zippedFastq%.*}"
	fileName="${zippedFastq%%.*}"
	
	# first unzip the file
	gunzip $zippedFastq
	
	# Map to hg19 with bowtie2
	echo "Mapping following file:"
	echo $fastqFile
	
	unsortedBAM="${fileName}_unsorted.bam"
	bowtie2 -x $hg19index -p 6 --time -U $fastqFile | samtools view -S -b - > $unsortedBAM
	
	#sort and index the bam file
	sorted="${fileName}_sorted"
	sortedBAM="${sorted}.bam"
	sortedBAI="${sortedBAM}.bai"
	samtools sort $unsortedBAM $sorted
	samtools index $sortedBAM
	
	# remove the unsorted BAM file
	rm $unsortedBAM

	#filter unmapped reads
	filteredBAM="${fileName}_filtered.bam"
	echo -e The numbers of mapped and unmapped reads, respectively:
	samtools view -c -F 4 $sortedBAM
	samtools view -c -f 4 $sortedBAM
	samtools view -h -F 4 -b $sortedBAM > $filteredBAM
	
	# remove unfiltered files
	rm $sortedBAM
	rm $sortedBAI

	#sort and index
	sorted="${fileName}_filtered_sorted"
	sortedBAM="${sorted}.bam"
	sortedBAI="${sortedBAM}.bai"
	samtools sort $filteredBAM $sorted
	samtools index $sortedBAM
	
	# remove unsorted file
	rm $filteredBAM

	#filter blacklisted sites & make index for resulting file
	finalBAM="${fileName}.bam"
	bedtools intersect -abam $sortedBAM -b $blacklistSites -v > $finalBAM
	samtools index $finalBAM
	
	# remove files without blacklist removed
	rm $sortedBAM
	rm $sortedBAI

	echo -e Final read count for this dataset after removal of blacklisted sequences:
	samtools view -c -F 4 $finalBAM
	
	# zip the fastq file back
	gzip $fastqFile
done







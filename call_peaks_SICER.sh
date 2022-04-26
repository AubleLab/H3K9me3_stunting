#!/bin/bash

# Provide a path to the SICER directory
echo -e "Provide a path to the directory that contains SICER.sh script."
read SICER

# get the input file
echo -e "Give a full name (including pathway) of an control BED file (coverage over input)."
read Input

# get the path to current directory
currentDirectory=pwd

# create output directory
mkdir SICER_output

# and get a path to that directory
outputDir="${currentDirectory}/SICER_output"


# Input variables: 
#["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] 
#["window size (bp)"] ["fragment size"] ["effective genome fraction"]   ["gap size (bp)"] ["FDR"]

# go through generated bed files (with coverage)
for bedFile in *.bed 
do 
	# call peaks
	sh $SICER/SICER.sh $currentDirectory $bedFile $Input $outputDir hg19 1 1000 425 0.74 5000 .01
	
done

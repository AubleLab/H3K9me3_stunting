#!/bin/bash
#the script takes the BAM file, normalizes it base on the normalization factor you give it and created bw file
# get chromosome sizes with refgenie
chromSize=$(refgenie seek hg19/fasta.chrom_sizes)

#the BAM file to be normalized
echo -e What is the name of BAM file?
read controlSample

#normalization factor : 1/(intended normalization factor: e.g. DESeq2 norm. fact)
echo -e What is the scale factor?
read X

#get everything before the first dot and use it as a file name
id="$( cut -d '.' -f 1 <<< "$controlSample" )"
suffix="_normalized.bg"
suffix_bw="_normalized.bw"

#get normalized bg file
bedtools genomecov -ibam $controlSample -bg -scale $X -g $chromSize > $id$suffix

#transform bg file to bw file
wigToBigWig -clip $id$suffix $chromSize $id$suffix_bw
#!/bin/bash
#the scrip goes through the BAM files 
#make directory for fastqc results
mkdir fastqc_results

#go through BAM files
for i in ./*.bam; do

id="$( cut -d '.' -f 2 <<< "$i" )"

echo ${id:2};
#run fastqc and put results in the above directory
fastqc -o ./fastqc_results/ $i

done

multiqc .


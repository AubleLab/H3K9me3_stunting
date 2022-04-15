# H3K9me3 in stunting
 All scripts relevant to the paper describing H3K27ac changes in stunted children.

## Data preprocessing
### 1) Map FASTQ files to hg19
###### Prerequisites:
+ [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
+ [refgenie](http://refgenie.databio.org/en/latest/)
+ [samtools](http://www.htslib.org/)
+ [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)

First check if you have bowtie2 index.\
`$ refgenie seek hg19/bowtie2_index`

If NOT found run following: \
`$ refgenie pull hg19/bowtie2_index` 

**From directory with compressed (.gz) FASTQ files run following script**\
`$ mapFASTQfiles_hg19.sh `

-Following sentence pops up:\
*`Give a full name (including pathway) of a file containing hg19 blacklisted sites.`*\
-Provide following (these are blacklisted sites defined by ENCODE - included here in */associated_files* folder):\
*`localPathTo_H3K27ac_in_stunting_folder/associated_files/hg19_blacklist.bed`*

### 2) Call peaks with SICER
###### Prerequisites:
+ [SICER](https://zanglab.github.io/SICER2/)

*********** edit **********
From directory with all of the hg19 BAM files run following:\
`$ callPeaks.sh`

In our settings we used input as a control for peak calling, when asked provide name of the input BAM file (including path name). !!! Input BAM file should be placed at different location from the other BAM files, otherwise peak-calling will be done also on this file).\
*`What is the name of the control dataset (including path)?`*\
e.g. *`localPathToInputBAM/input.bam`*\
\
Outputs from peak-calling are placed into individual folders named after individual BAM files.

*********** edit **********


### 3) Create count tables
Move the output files from previous step ending with "_Enhancer.bed" into a separate folder and from here run following (separate for each age group): \
`cat *.bed | sort -k1,1 -k2,2n | bedtools merge -i stdin > masterPeaks.bed`

Place the masterEnhancer.bed file into a folder containing all BAM files for a given age group and from this folder run: \
`bedtools multicov -bams *.bam -bed masterPeaks.bed > countTable.txt`

First three columns in the table are genomic coordinates of regions of interest followed by counts for individual samples. Make sure to add sample names to the sample columns based on their order within the folder (`ls *.bam`).

### 4) Map FASTQ files to dm6
First check if you have bowtie2 index.\
`$ refgenie seek dm6/bowtie2_index`

If NOT found run following: \
`$ refgenie pull dm6/bowtie2_index` 

**From directory with compressed (.gz) FASTQ files run following script**\
`$ mapFASTQfiles_dm6.sh `

-Following sentence pops up:\
*`Give a full name (including pathway) of a file containing dm6 blacklisted sites.`*\
-Provide following (these are blacklisted sites defined by ENCODE - included here in */associated_files* folder):\
*`localPathTo_H3K27ac_in_stunting_folder/associated_files/dm6_blacklist.bed`*

### 4) Get heatmaps and average profile over genes and peaks
###### Prerequisites:
+ [deepTools](https://deeptools.readthedocs.io/en/develop/)

*********** edit number **********
### XXX) Create database for LOLA from CISTROME database
The CISTROME database for human transcription factors and histone marks can be downloaded from [CISTROME website](http://cistrome.org/db/#/bdown).\
From within a folder containing info about downloaded files: \
*human_factor_full_QC.txt*, \
*human_hm_full_QC.txt*,\
and actual folders with downloaded BED files: \
*human_factor/*, \
*human_hm/* 

Run following R script:\
`filterBloodSpecificFiles.R`\
\
The script creates new directories, where only BED files coming from blood cells are placed (*human_factor_blood/*, *human_hm_blood/*) and new info sheets are created (*human_factor_blood_QC.txt*, *human_hm_blood_QC.txt*).

LOLA takes database input in form of GRangesList object. To convert BED files into GRangesList objects, run following script from within directory containing the *human_factor_blood/*, and *human_hm_blood/* directories.\
`makeGRangesList_forLOLA.R`\
\
The script creates 2 GRangesList objects from the BED files: *human_factor_blood.Rdata*, and *human_hm_blood.Rdata*. It also produces a list of empty BED files (*human_factor_blood_unloaded.txt*, *human_hm_blood_unloaded.txt*) that were not added to the GRangesList objects, and will be therefore excluded from an annotation table created in the following step. 

LOLA annotation file requires number of lines within each bed file. From terminal run following: \
`$ cd localPathTo/human_factor_blood` \
`$ num_of_lines_in_bed_files.sh > human_factor_blood_bedSize.txt` \
`$ mv human_factor_blood_bedSize.txt ..`

`$ cd localPathTo/human_hm_blood` \
`$ num_of_lines_in_bed_files.sh > human_hm_blood_bedSize.txt` \
`$ mv human_hm_blood_bedSize.txt ..`

The region and collection annotations for LOLA can then be created by running following R script from within directory containing the *human_factor_blood/*, and *human_hm_blood/* directories.
`makeAnnotFiles_LOLA.R` 

The resulting files are then *human_factor_blood_regionAnno.csv*, *human_factor_blood_collectionAnno.csv*, *human_hm_blood_regionAnno.csv*, and *human_hm_blood_collectionAnno.csv*. 

You can now move the following list of files to a separate directory, where you want to have LOLA database stored, and you can remove the rest of the files.
Files for LOLA database: \
*human_factor_blood_collectionAnno.csv,\
human_factor_blood_regionAnno.csv,\
human_factor_blood.Rdata,\
human_hm_blood_collectionAnno.csv,\
human_hm_blood_regionAnno.csv,\
human_hm_blood.Rdata*

You can then run LOLA with help of following R script, where information about location of created database, location of BED file of interest and location of universe BED file must be provided between lines 9-25. For more information on use od LOLA, you can click on the next [link](http://databio.org/lola/).

`LOLA_cistrome.R`


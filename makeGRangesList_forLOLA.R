# Script takes bed files from a folder - make GRangesList out of them
# when there is an empty bed file - the script skips it and writes down the 
# bed file, so it doesn't appear in the annotation table after that

rm(list = ls())
library(GenomicRanges)
library(tidyverse)

for (folder in c("human_hm_blood", "human_factor_blood")){
  # go through the files within the folder - 
  # save them into GRangesList
  bedFiles = list.files(folder)
  bedFilesFull = list.files(folder, full.names = T)

  # set up variable for empty bed files
  faultyBed = c()
  for (i in bedFilesFull){
    
    # skip empty bed files, but write down the record
    if (file.info(i)$size == 0){
      faultyBed = c(faultyBed, i)
      next
    }
    
    bedFile = read.delim(i, header = F) %>% 
      select(V1:V3)
    colnames(bedFile) = c("chr", "start", "end")
    bedFile = makeGRangesFromDataFrame(bedFile)
    bedFile = keepStandardChromosomes(bedFile, pruning.mode = "coarse")
    bedFile = GRangesList(bedFile)
    
    if (i == bedFilesFull[1]){
      regionGRL = bedFile
    } else {
      regionGRL = c(regionGRL, bedFile)
    }
    rm(bedFile)
  }
  # save the GRangesList with name corresponding to the originating folder
  save(regionGRL, file = paste0(folder, ".Rdata"))
  # save list of empty BED files - these are not in the GRangesList
  write.table(faultyBed, paste0(folder, "_unloaded.txt"), quote = F, row.names = F, col.names = F, sep = "\n")
  
}



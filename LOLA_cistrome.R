rm(list = ls())

library(LOLA)
library(GenomicRanges)
library(tidyverse)
library(simpleCache)
library(data.table)

##### PROVIDE following ######
# give path, where results should be stored
resultPath = "path/to/results"

# give path to the directory with created files from downloaded CISTROME database
dbLocation = "path/to/CISTROME_for_LOLA/database"

# select which database to use - histone mark or transcription factor (just comment the other)
element = "human_factor_blood"
#element = "human_hm_blood"

# provide path to universe (including the file name) - file must be hg38!!!
pathToUniverse = "path/to/bed/file/with/full/set/of/regions"

# provide path to set (including the file name) = regions of interest - file must be hg38!!!
pathToSet = "path/to/bed/file/with/set/of/regions/of/interest"
##################################

load(paste0(dbLocation, "/", element, ".Rdata"))
fileAnnotation = read.csv(paste0(dbLocation, "/", element, "_regionAnno.csv"))
regionAnno = fileAnnotation %>% 
  unite(description, c(GSMID, Cell_line), sep = "_") %>% 
  mutate(cellType = Cell_type) %>% 
  mutate(tissue = Tissue_type) %>% 
  mutate(dataSource = NA) %>% 
  mutate(antibody = Factor) %>% 
  mutate(treatment = NA) %>% 
  select(filename, cellType, description, tissue, dataSource, antibody, treatment, collection, size)
regionAnno = as.data.table(regionAnno)

collectionAnno = fread(paste0(dbLocation, "/", element, "_collectionAnno.csv"))

regionDB= list(dbLocation = dbLocation, 
                   regionAnno = regionAnno, 
                   collectionAnno = collectionAnno, 
                   regionGRL= regionGRL)

# upload peaks for universe
bed = read.delim(pathToUniverse, header = F) %>% 
  dplyr::rename(chr = V1) %>% 
  dplyr::rename(start = V2) %>% 
  dplyr::rename(end = V3) %>% 
  mutate(strand = "*")

# transform to GRanges
universe = makeGRangesFromDataFrame(bed)
universe = keepStandardChromosomes(universe, pruning.mode = "coarse")

# get significant peaks
set = read.delim(pathToSet, header = F) %>% 
  dplyr::rename(chr = V1) %>% 
  dplyr::rename(start = V2) %>% 
  dplyr::rename(end = V3) %>% 
  mutate(strand = "*")

set = makeGRangesFromDataFrame(set)
set = keepStandardChromosomes(set, pruning.mode = "coarse")


#run LOLA
locResults = runLOLA(set, universe, regionDB, cores=1)
write.csv(locResults, paste0(resultPath,element, "_allResults.csv"), row.names = F)


# filter out significant results
sigLOLA = locResults %>% 
  filter(qValue <= 0.05)
write.csv(sigLOLA, paste0(resultPath, element, "_significantResults.csv"), row.names = F)







rm(list = ls())
library(tidyverse)

# 1) human_factor_blood

type = "human_factor_blood"

# remove the empty BED files from list of BED files
bedFiles = list.files(type)
emptyBedFiles = read.delim(paste0(type, "_unloaded.txt"), header = F, sep = "\n")
emptyBedFiles = basename(emptyBedFiles$V1)
usedBedFiles = bedFiles[!(bedFiles %in% emptyBedFiles)]

# bed file sizes: 
sizes = read.delim(paste0(type, "_bedSize.txt"), header = F)
filename = sizes[seq(1, nrow(sizes), by = 2),]
size = sizes[seq(2, nrow(sizes), by = 2),]
sizeTable = data.frame(filename, size)

# merge cistrome annotation with used bed file names 
cistromeAnnotation = read.delim(paste0(type, "_QC.txt")) %>% 
  mutate(DCid = as.character(DCid))

usedBedFilesTable = data.frame(filename = usedBedFiles) %>% 
  separate(filename, c("DCid", "toss"), sep = "_", remove = F, extra = "merge") %>% 
  select(-toss)

regionAnno = usedBedFilesTable %>% 
  left_join(cistromeAnnotation) %>% 
  left_join(sizeTable) %>% 
  mutate(collection = type)

# write the annotation into a CSV
write.csv(regionAnno, paste0(type, "_regionAnno.csv"), row.names = F)

# create collection annotaion:
collectionAnno = data.frame(collectionname = type, 
                            collector = "YourName", 
                            date = Sys.Date(), 
                            source = "cistrome_factor", 
                            description = "http://cistrome.org/db/#/bdown")
write.csv(collectionAnno, "human_factor_blood_collectionAnno.csv", quote = F, row.names = F)


# 2) human_hm_blood

type = "human_hm_blood"

# remove the empty BED files from list of BED files
bedFiles = list.files(type)
emptyBedFiles = read.delim(paste0(type, "_unloaded.txt"), header = F, sep = "\n")
emptyBedFiles = basename(emptyBedFiles$V1)
usedBedFiles = bedFiles[!(bedFiles %in% emptyBedFiles)]

# bed file sizes: 
sizes = read.delim(paste0(type, "_bedSize.txt"), header = F)
filename = sizes[seq(1, nrow(sizes), by = 2),]
size = sizes[seq(2, nrow(sizes), by = 2),]
sizeTable = data.frame(filename, size)

# merge cistrome annotation with used bed file names 
cistromeAnnotation = read.delim(paste0(type, "_QC.txt")) %>% 
  mutate(DCid = as.character(DCid))

usedBedFilesTable = data.frame(filename = usedBedFiles) %>% 
  separate(filename, c("DCid", "toss"), sep = "_", remove = F, extra = "merge") %>% 
  select(-toss)

regionAnno = usedBedFilesTable %>% 
  left_join(cistromeAnnotation) %>% 
  left_join(sizeTable) %>% 
  mutate(collection = type)
# write the annotation into a CSV
write.csv(regionAnno, paste0(type, "_regionAnno.csv"), row.names = F)

# create collection annotaion:
collectionAnno = data.frame(collectionname = type, 
                            collector = "YourName", 
                            date = Sys.Date(), 
                            source = "cistrome_hm", 
                            description = "http://cistrome.org/db/#/bdown")
write.csv(collectionAnno, "human_hm_blood_collectionAnno.csv", quote = F, row.names = F)



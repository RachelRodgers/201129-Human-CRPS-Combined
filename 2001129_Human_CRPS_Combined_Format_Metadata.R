# 201129_Human_CRPS_Combined_Format_Metadata.R

# For this analysis, we are combining human CRPS data from across five sequencing runs.
#   Here I will be combining the information from the mapping files.

library("tidyverse")

mappingFileList <- c("run1" = "mappingFile_190201.txt",
                     "run2" = "mappingFile_190318.txt",
                     "run3" = "mappingFile_190523.txt",
                     "run4" = "mappingFile_200127.txt",
                     "run5" = "mappingFile_201019_LaraC_V4_Mouse_Human_CRPS.txt")

mappingFiles <- map(paste0("../documents/", mappingFileList), read.delim, 
                    header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                    check.names = FALSE)

names(mappingFiles) <- names(mappingFileList)


# Find common column names among files
allColNames <- map(mappingFiles, names)

commonColNames <- Reduce(f = intersect, x = allColNames)

# For each mapping file, pull out these common columns and add a run column
trimmedMapFileList <- vector(mode = "list", length = length(mappingFiles))
  
for (i in 1:length(mappingFiles)) {
  
  currentRun <- names(mappingFiles)[i]
  currentMap <- mappingFiles[[i]]
  
  trimmedDF <- currentMap %>% 
    select(all_of(commonColNames)) %>% 
    mutate("run" = currentRun) %>% 
    select(`#SampleID`, BarcodeSequence, LinkerPrimerSequence, run, Description)
  
  trimmedMapFileList[[i]] <- trimmedDF
  names(trimmedMapFileList)[i] <- currentRun
  
}

# Combine the trimmed mapping files
combinedMap <- Reduce(f = rbind, x = trimmedMapFileList)

# Remove "CRPS Mouse", format the Description column, add other needed columns
table(combinedMap$Description)

combinedMapFmt <- combinedMap %>% 
  filter(Description != "CRPS mouse") %>% 
  mutate(Description = trimws(Description),
         Description = str_replace(Description, 
                                   pattern = " ", replacement = "_"),
         Sample_Type = case_when(Description %in% c("Acute", "Chronic") ~ "CRPS",
                                 Description == "Control" ~ "HHC",
                                 TRUE ~ Description),
         Sample_Status = case_when(Description %in% c("Acute", "Chronic") ~ "CRPS",
                                   Description %in% c("Control", "BioBank_Control") ~ "healthy",
                                   TRUE ~ Description)) %>% 
  select(`#SampleID`, BarcodeSequence, LinkerPrimerSequence, run, Sample_Type,
         Sample_Status, Description)

table(combinedMapFmt$Description)

write.table(combinedMapFmt, file = "../data/mappingFile_201129_Human_CRPS_Combined.txt",
            quote = FALSE, sep = "\t", row.names = FALSE)

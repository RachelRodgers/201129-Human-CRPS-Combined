# 201129_Human_CRPS_Combined_SampleDataCleaning.R

# Cleaning up the sample data that will be used for the 16S analysis.

library("tidyverse")

sampleData <- read.delim("../data/sampleDataLC.txt", stringsAsFactors = FALSE,
                         check.names = FALSE)

# remove columns:
sampleDataModified <- sampleData %>% 
  select(-c("BarcodeSequence", "LinkerPrimerSequence", "HTCF_path",
            "index_1", "index_2", "bc_1", "bc_2", contains("16S.Barcode")))

# change case of column names and replace "." with "_":
colNames <- names(sampleDataModified)
colNamesLower <- colNames %>% 
  str_to_lower() %>% 
  str_replace(pattern = "\\.", replacement = "_")

names(sampleDataModified) <- colNamesLower

# add a "household" column to more easily tell which samples are matched
crpsHousehold <- sampleDataModified %>% 
  filter(sample_type == "CRPS") %>% 
  select(sample) %>% 
  mutate("number" = str_extract(string = sample, pattern = "[:digit:]+"),
         "letter" = str_extract(string = sample, pattern = "[:alpha:]?[:digit:]?$"),
         "household" = paste(number, letter, sep = "-"),
         "household" = str_remove(string = household, pattern = "[:digit:]?$")) %>% 
  select(sample, household) %>% 
  deframe()

hhcHousehold <- sampleDataModified %>% 
  filter(sample_type == "HHC") %>% 
  select(sample, sample_match) %>% 
  mutate("number" = str_extract(string = sample_match, pattern = "[:digit:]+"),
         "letter" = str_extract(string = sample_match, pattern = "[:alpha:]?[:digit:]?$"),
         "household" = paste(number, letter, sep = "-"),
         "household" = str_remove(string = household, pattern = "[:digit:]?$")) %>% 
  select(sample, household) %>% 
  deframe()

householdLUT <- c(crpsHousehold, hhcHousehold)

sampleDataModified <- sampleDataModified %>% 
  mutate("household" = ifelse(sample %in% names(householdLUT),
                              yes = householdLUT[sample],
                              no = sample_match)) %>% 
  select("sample", "sample_type", "sample_status", "description", "run",
         "condition", "crps_match", "sample_match", "household", everything())

saveRDS(sampleDataModified, "../data/sampleDataModified.RDS")

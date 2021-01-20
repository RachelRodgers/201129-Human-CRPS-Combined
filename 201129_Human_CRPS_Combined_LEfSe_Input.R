# 201129_Human_CRPS_Combined_LEfSe_Input.R

# Generate LEfSe Input Data

# Comparisons to be conducted:
#   1.) Acute vs Chronic CRPS
#   2.) Acute CRPS vs Acute BioBank
#   3.) Chronic CRPS vs Chronic BioBank
#   4.) Acute CRPS vs Acute HHC
#   5.) Chronic CRPS vs Chronic HHC
#   6.) All CRPS vs All HHC
#   7.) All CRPS vs All BioBank

source("./shared_R_scripts/GenerateLefseData.R")
source("./shared_R_scripts/Helper_Functions.R")

physeqBacteria <- readRDS("../data/physeqObjects/physeqBacteria.RDS")

sampleData <- as(sample_data(physeqBacteria), "data.frame")

conditionNames <- c("crps_acute_chronic",
                    "acute_crps_biobank",
                    "chronic_crps_biobank",
                    "acute_crps_hhc",
                    "chronic_crps_hhc",
                    "crps_hhc",
                    "crps_biobank")

# What column to find the condition names in:
columnVars <- c(rep("crps_match", 5), rep("sample_type", 2))
names(columnVars) <- conditionNames
columnVars <- as.list(columnVars)

# What labels meet each condition:
conditionVecs <- list(c("CRPSA", "CRPSC"),
                      c("CRPSA", "NRAM"),
                      c("CRPSC", "NRCM"),
                      c("CRPSA", "HHCA"),
                      c("CRPSC", "HHCC"),
                      c("CRPS", "HHC"),
                      c("CRPS", "BioBank_Control"))
names(conditionVecs) <- conditionNames

# initialize empty list to hold lefse input data
lefseList <- vector(mode = "list", length = length(conditionNames))

for (i in 1:length(conditionNames)) {

  currCondition <- conditionNames[i]
  currColumnVar <- columnVars[[currCondition]]
  currConditionVec <- conditionVecs[[currCondition]]

  keepSamples <- sampleData %>% 
    filter(!!(sym(currColumnVar)) %in% currConditionVec) %>% 
    pull(sample)
  
  psSubset <- prune_samples(keepSamples, physeqBacteria) %>% 
    RemoveMissingTaxa()
    
  lefse <- GenerateLefseData(physeqObject = psSubset,
                             categoryColumnName = currColumnVar,
                             sampleColumnName = "sample")
  
  lefseList[[i]] <- lefse
  names(lefseList)[i] <- currCondition
}

# save data
paths <- paste0("../data/LEfSe_input/", names(lefseList), ".txt")

walk2(.x = lefseList, .y = paths,
      .f = ~ write.table(x = .x, file = .y, quote = FALSE,
                         sep = "\t", row.names = TRUE, col.names = FALSE))

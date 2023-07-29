# 201129_HumanCRPSCombined_RDPandSILVA_FormatPhyloseq.R

# Format the RDP and SILVA phyloseq objects needed for a side-by-side analysis.

library("phyloseq")
library("tidyverse")

#----- Generate Combined Metadata from RDP & SILVA -----#

# sample data was originally formatted separately for each data set
sampleDataRDP <- readRDS("../data/RDP/RDataObjects/sampleDataFinalRDP.RDS")
sampleDataSILVA <- readRDS("../data/SILVA/RDataObjects/sampleDataFinalSILVA.RDS")

# select columns to combine together
silvaCols <- colnames(sampleDataSILVA)
rdpCols <- colnames(sampleDataRDP)

# prevent dupliated columns:
# in the combined sample data, only keep columns from SILVA not already in RDP
sampleDataCombined <- sampleDataSILVA %>% 
  select(sample, silvaCols[!silvaCols %in% rdpCols]) %>% 
  merge(sampleDataRDP, by = "sample", all = FALSE) %>%
  mutate(bad_household_chr = ifelse(bad_household == TRUE,
                                    yes = "Bad Household",
                                    no = "Good Household"))

#----- RDP v16-----#
# originally 69 samples
physeqRawRDP <- readRDS("../data/RDP/RDataObjects/physeqObjects/ps0.rdp_single.RDS")
sampleDataRDP <- sampleDataCombined
row.names(sampleDataRDP) <- sampleDataRDP$sample
physeqMergedRDP <- merge_phyloseq(physeqRawRDP, sample_data(sampleDataRDP))
physeqMergedRDP # 2115 x 67

rdpMissing <- sample_names(physeqRawRDP)[!sample_names(physeqRawRDP) %in% sample_names(physeqMergedRDP)]
# sample Baldridge007B and Baldridge Control have been correctly dropped

#----- SILVA v138.1 -----#
# originally 118 samples
physeqRawSILVA <- readRDS("../data/SILVA/RDataObjects/physeqObjects/ps0.silva_single.RDS")
sampleDataSILVA <- sampleDataCombined
row.names(sampleDataSILVA) <- sampleDataCombined$sequencing_file_name
physeqMergedSILVA <- merge_phyloseq(physeqRawSILVA,
                                    sample_data(sampleDataSILVA))
physeqMergedSILVA # 2668 x 67 samples

# These phyloseq objects are ready for analysis:

# 1.) RDP: 2115 x 67, use "sample" as the sample name
saveRDS(physeqMergedRDP,
        "../data/RDP/RDataObjects/physeqObjects/physeqMergedRDP.RDS")

# 2.) SILVA: 2668 x 67, use "sequencing_file_name" as the sample name
saveRDS(physeqMergedSILVA,
        "../data/SILVA/RDataObjects/physeqObjects/physeqMergedSILVA.RDS")

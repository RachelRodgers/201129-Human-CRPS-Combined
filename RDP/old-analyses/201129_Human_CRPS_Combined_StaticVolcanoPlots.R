# 201129_Human_CRPS_Combined_StaticVolcanoPlots.R

library("phyloseq")
library("gridExtra")
library("ggrepel")
library("ggpubr")
library("DESeq2")
library("tidyverse")

load("201129_Human_CRPS_Combined_16S_Analysis.RData")
source("./shared_R_scripts/Helper_Functions.R")

#----- Household Controls vs. CRPS -----#

physeqHousehold <- physeqSubsetList$household
View(as(sample_data(physeqHousehold), "matrix"))

# Run DESeq analysis on Household
ddsHousehold <- phyloseq_to_deseq2(physeqHousehold, 
                                   design = as.formula("~ sample_type"))
countsHousehold <- counts(ddsHousehold)
geoMeans <- apply(countsHousehold, 1, function(row) {
  if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
})     

ddsHousehold <- estimateSizeFactors(ddsHousehold, geoMeans = geoMeans)
ddsAnalysisHousehold <- DESeq(ddsHousehold, test = "Wald", fitType = "local",
                              betaPrior = FALSE)
ddsResultsHousehold <- results(ddsAnalysisHousehold,
                               contrast = c("sample_type",
                                            "CRPS", "HHC"))
mcols(ddsResultsHousehold)

# household results table
resTableHousehold <- GenerateDESeqResultsTable(physeqHousehold,
                                               ddsResultsHousehold)

householdVolcano <- PlotStaticVolcano(physeqHousehold, resTableHousehold, 
                                      sigThreshold = 0.05)

#----- Chronic vs. Acute CRPS -----#

View(as(sample_data(physeqBacteria), "matrix"))

physeqCRPS <- physeqBacteria %>% 
  subset_samples(sample_status == "CRPS") %>% 
  RemoveMissingTaxa()
physeqCRPS

# run DESeq2 analysis on CRPS
ddsCRPS <- phyloseq_to_deseq2(physeqCRPS, design = as.formula("~ description"))
countsCRPS <- counts(ddsCRPS)

geoMeansCRPS <- apply(countsCRPS, 1, function(row) {
  if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
})

ddsCRPS <- estimateSizeFactors(ddsCRPS, geoMeans = geoMeansCRPS)

ddsAnalysisCRPS <- DESeq(ddsCRPS, test = "Wald", fitType = "local",
                         betaPrior = FALSE)

#ddsResultsCRPS <- results(ddsAnalysisCRPS,
                          #contrast = c("description", "Acute", "Chronic"))
ddsResultsCRPS <- results(ddsAnalysisCRPS)
mcols(ddsResultsCRPS)

# CRPS results table
resTableCRPS <- GenerateDESeqResultsTable(physeqCRPS, ddsResultsCRPS)

crpsVolcano <- PlotStaticVolcano(physeqCRPS, resTableCRPS, sigThreshold = 0.05)
crpsVolcano

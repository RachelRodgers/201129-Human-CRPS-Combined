# 201129_Human_CRPS_Combined_StaticBetaPlots.R

library("phyloseq")
library("vegan")
library("gridExtra")
library("knitr")
library("plotly")
library("microbiome")
library("data.table")
library("ggrepel")
library("ggpubr")
library("scales")
library("plyr")
library("DT")
library("picante")
library("gridExtra")
library("tidyverse")

load("201129_Human_CRPS_Combined_16S_Analysis.RData")
source("./shared_R_scripts/Helper_Functions.R")

#----- CRPS vs. BioBank (BioBank) -----#

physeqBioBank <- physeqBacteria %>% 
  subset_samples(sample_type %in% c("BioBank_Control", "CRPS")) %>% 
  RemoveMissingTaxa()

set.seed(13905870)
ordBioBank <- pmap(.l = list(distance = list("unweighted" = "unifrac", 
                                             "weighted" = "wunifrac")),
                   .f = ordinate, 
                   physeq = physeqBioBank, method = "PCoA")

PlotBioBank <- function(currentOrd) {
  plot_ordination(physeqBioBank,
                  ordination = currentOrd,
                  color = "sample_type", shape = "sample_type") +
    scale_shape_manual(values = c(16, 15)) +
    scale_color_manual(values = c("#FECC66", "#FB0106")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(size = 4) +
    stat_ellipse(type = "norm")
}

pmap(.l = list(currentOrd = ordBioBank), .f = PlotBioBank)

#----- Acute vs. Chronic CRPS (CRPS) -----#

physeqCRPS <- physeqBacteria %>% 
  subset_samples(description %in% c("Acute", "Chronic")) %>% 
  RemoveMissingTaxa()

set.seed(13905870)
ordCRPS <- pmap(.l = list(distance = list("unweighted" = "unifrac",
                                          "weighted" = "wunifrac")),
                   .f = ordinate, 
                   physeq = physeqCRPS, method = "PCoA")

PlotCRPS <- function(currentOrd) {
  plot_ordination(physeqCRPS,
                  ordination = currentOrd,
                  color = "description", shape = "description") +
    scale_shape_manual(values = c(17, 18)) +
    scale_color_manual(values = c("#CC66FF", "#408002")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(size = 4) +
    stat_ellipse(type = "norm")
}

pmap(.l = list(currentOrd = ordCRPS), .f = PlotCRPS)

#----- CRPS vs HHC (household) -----#

set.seed(13905870)
ordHH <- pmap(.l = list(distance = list("unweighted" = "unifrac",
                                        "weighted" = "wunifrac")),
                .f = ordinate, 
                physeq = physeqSubsetList$household, method = "PCoA")

PlotHH <- function(currentOrd) {
  plot_ordination(physeqSubsetList$household,
                  ordination = currentOrd,
                  color = "sample_type", shape = "sample_type") +
    scale_shape_manual(values = c(15, 19)) +
    scale_color_manual(values = c("#FB0106", "#3333CC")) +
    theme_bw() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_point(size = 4) +
    stat_ellipse(type = "norm")
}

pmap(.l = list(currentOrd = ordHH), .f = PlotHH)

set.seed(13905870)
RunAdonis(physeqSubsetList$household, "sample_type", "wunifrac")


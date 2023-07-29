# 201129-Human-CRPS-Combined

This repository contains the code and associated data for the analysis of the human CRPS data using both RDP and SILVA annotated data sets. Please see the resources directory for input and intermediate data, as well as ready-to-analyze phyloseq objects (for use with the script **201129_HumanCRPSCombined_RDP-SILVA-16S-Analysis.Rmd**. To rerun the analysis pipeline, please run the scripts in the order listed below. 

---

## Getting Started:

Run the following command to download all code, including submodules, and associated data:

```
git clone --recurse-submodules git@github.com:RachelRodgers/201129-Human-CRPS-Combined.git
```
Please note, data and documents subdirectories are referenced as being in separate data/ and documents/ directories on the same level as this repository (as opposed to a resources/ subdirectory as it appears when cloned). This is because the total amount of data associated with this project is not sourced controlled on GitHub, and the minimal required data has been included here as a convenience. Please move the data and documents directories to the same level as this repo, or edit the relative paths accordingly.

## RDP-Specific Scripts:
1. 201129_HumanCRPSCombined_RDP-Format-Demuxing-Metadata.Rmd
2. 201129_HumanCRPSCombined_RDP-DADA2.Rmd
3. 201129_HumanCRPSCombined_RDP-Format-16S-Analysis-Metadata.Rmd

## SILVA-Specific Scripts:
1. 201129_HumanCRPSCombined_SILVA-Format-Metadata.Rmd
2. 201129_HumanCRPSCombined_SILVA-DADA2.Rmd

## Combined Scripts:
1. 201129_HumanCRPSCombined_RDP-SILVA-Format-Phyloseq.R
2. 201129_HumanCRPSCombined_RDP-SILVA-16S-Analysis.Rmd

---

## Getting Help
Please open an issue or contact rachel.rodgers@wustl.edu.

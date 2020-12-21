# iCPAGdb

iCPAGdb is designed to facilitate rapid analysis of genetic correlation across thousands of GWAS simultaneously. Current verion of iCPAGdb contains GWAS summary statistic from >4400 diseases/traits, and allows users to explore pre-computed correlations across all existing diseases and/or upload their own GWAS to identify and explore shared SNPs between their own GWAS and >4400 diseases/traits.

This repo contains two parts
1) python code for iCPAGdb

2) R shiny codes for Web browser
Resources for the iCPAGdb web app

# Quick start

## download PLINK and database
1) download PLINK 1.9 () and place to folder "plink_bins"
2) download database file (~15Gb) to "db" folder from Dropbox (https://www.dropbox.com/sh/na23jflxcgk0nib/AAAoOj3gB8k8j_dNH1UBFzeZa?dl=0).

The following folder structure contains all required codes and data file:
pyCPAGdb 
├── _utils.py
├── anno_parent_efo.py
├── main.py
├── stats.py
├── plink_bins
│   ├── plink
│   └── prettify
├── db
│   ├── cpag_gwasumstat_v1.1.db
│   ├── cpag_gwasumstat_v1.1_AFR.ld0.4.db
│   ├── cpag_gwasumstat_v1.1_EAS.ld0.4.db
│   ├── cpag_gwasumstat_v1.1_EUR.ld0.4.db
│   ├── gwas-efo-trait-mappings.txt
│   └── lddat
│       ├── AFR_1kg_20130502_maf01.bed
│       ├── AFR_1kg_20130502_maf01.bim
│       ├── AFR_1kg_20130502_maf01.fam
│       ├── EAS_1kg_20130502_maf01.bed
│       ├── EAS_1kg_20130502_maf01.bim
│       ├── EAS_1kg_20130502_maf01.fam
│       ├── EUR_1kg_20130502_maf01.bed
│       ├── EUR_1kg_20130502_maf01.bim
│       └── EUR_1kg_20130502_maf01.fam



## run example


iCPAGdb 

Web browser: http://cpag.oit.duke.edu/ 

This repo contains all codes of the webbrowers.

iCPAGdb is designed to facilitate rapid analysis of genetic correlation across thousands of GWAS simultaneously. Current verion of iCPAGdb contains GWAS summary statistic from >4400 diseases/traits, and allows users to explore pre-computed correlations across all existing diseases and/or upload their own GWAS to identify and explore shared SNPs between their own GWAS and >4400 diseases/traits.

This repo contains two parts
1) python code for iCPAGdb

2) R shiny codes for Web browser
Resources for the iCPAGdb web app

# Quick start

## Download PLINK and database
1) download PLINK 1.9 () and place to folder "plink_bins"  <br/>
2) download database file (~15Gb) to "db" folder from Dropbox (https://www.dropbox.com/sh/na23jflxcgk0nib/AAAoOj3gB8k8j_dNH1UBFzeZa?dl=0).

The following folder structure contains all required codes and data file:<br/>
"pyCPAGdb <br/>
├── _utils.py <br/>
├── anno_parent_efo.py <br/>
├── main.py <br/>
├── stats.py <br/>
├── plink_bins <br/>
│   ├── plink <br/>
│   └── prettify <br/>
├── db <br/>
│   ├── cpag_gwasumstat_v1.1.db <br/>
│   ├── cpag_gwasumstat_v1.1_AFR.ld0.4.db <br/>
│   ├── cpag_gwasumstat_v1.1_EAS.ld0.4.db <br/>
│   ├── cpag_gwasumstat_v1.1_EUR.ld0.4.db <br/>
│   ├── gwas-efo-trait-mappings.txt <br/>
│   └── lddat <br/>
│       ├── AFR_1kg_20130502_maf01.bed <br/>
│       ├── AFR_1kg_20130502_maf01.bim <br/>
│       ├── AFR_1kg_20130502_maf01.fam <br/>
│       ├── EAS_1kg_20130502_maf01.bed <br/>
│       ├── EAS_1kg_20130502_maf01.bim <br/>
│       ├── EAS_1kg_20130502_maf01.fam <br/>
│       ├── EUR_1kg_20130502_maf01.bed <br/>
│       ├── EUR_1kg_20130502_maf01.bim <br/>
│       └── EUR_1kg_20130502_maf01.fam <br/>

## Run example

 <br/>

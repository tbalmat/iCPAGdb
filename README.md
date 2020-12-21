# iCPAGdb 

Web browser: http://cpag.oit.duke.edu/ 

This repo contains all codes of the web browser.

iCPAGdb is designed to facilitate rapid analysis of genetic correlation across thousands of GWAS simultaneously. Current verion of iCPAGdb contains GWAS summary statistic from >4400 diseases/traits, and allows users to explore pre-computed correlations across all existing diseases and/or upload their own GWAS to identify and explore shared SNPs between their own GWAS and >4400 diseases/traits.

This repo contains two parts
1) python3 code for iCPAGdb.

2) R shiny code for Web browser


# Quick start

## Download PLINK and database
1) download PLINK 1.9 () and place to folder "plink_bins"  <br/>
2) download database file (~15Gb) to "db" folder from Dropbox (https://www.dropbox.com/sh/na23jflxcgk0nib/AAAoOj3gB8k8j_dNH1UBFzeZa?dl=0).

The following folder structure contains all required codes and data file:<br/>
pyCPAGdb <br/>
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

## Run example (Shell command)

### example 1

Serum metabolites/xenobiotics (Shin et al. 2014) vs. Human disease 

```python3 main.py cpagdb --threads 2 --subtype NHGRI --NHGRI-Pcut 5e-8 --subtype BloodMetabolites,BloodXenobiotic --Pcut 1e-5 --lddb-pop EUR --outfile NHGRI-p1e-05-BloodMetabolitesXenobiotic-p1e-05-EUR.csv```

annotate phenotype:

```python3 main.py post_analysis --anno-ontology --anno-cols Trait1 --infile output/NHGRI-p1e-05-BloodMetabolitesXenobiotic-p1e-05-EUR.csv --outfile NHGRI-p1e-05-BloodMetabolitesXenobiotic-p1e-05-EUR.csv```

### example 2
```python3 main.py cpagdb --threads 2 --subtype H2P2 --H2P2-Pcut 1e-7 --lddb-pop EUR --outfile output/H2P2-p1e-07-EUR.csv```

### example 3 (user GWAS)

download example file (iCPAGdb-Sample-GWAS-top_EllinghausPCs_covid19.csv) from "Upload and compute CPAG" page at [here](http://cpag.oit.duke.edu/)

```python main.py usr-gwas --threads 10 --infile iCPAGdb-Sample-GWAS-top_EllinghausPCs_covid19.csv --SNPcol "avsnp150" --delimitor "," --Pcol "p_value" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --ld-clump 1 --outfile top_EllinghausPCs_covid19_pcut1e-5_icpagdb_out.csv```

annotate phenotype:

```python3 main.py post_analysis --anno-ontology --anno-cols Trait2 --infile top_EllinghausPCs_covid19_pcut1e-5_icpagdb_out.csv --outfile top_EllinghausPCs_covid19_pcut1e-5_icpagdb_out_addEFO.csv```



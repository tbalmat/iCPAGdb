
#python main.py cpagdb --threads 30 --subtype H2P2 --H2P2-Pcut 1e-3
#python main.py cpagdb --threads 30 --subtype H2P2 --H2P2-Pcut 1e-5
#python main.py cpagdb --threads 30 --subtype H2P2 --H2P2-Pcut 1e-6
#python main.py cpagdb --threads 30 --subtype H2P2 --H2P2-Pcut 1e-7
#python main.py cpagdb --threads 10 --subtype H2P2 --H2P2-Pcut 1e-4

#python main.py cpagdb --threads 6 --subtype mol_gwas --H2P2-Pcut 1e-5 --Pcut 5e-8 --outfile "./output/cpag_output_molGWASpcut5e-8_(h2p2pcut1e-5)_out.csv"
#python main.py cpagdb --threads 6 --subtype mol_gwas --H2P2-Pcut 1e-4 --Pcut 5e-8 --outfile "./output/cpag_output__molGWASpcut5e-8_(h2p2pcut1e-4)_out.csv"
#python main.py cpagdb --threads 6 --subtype mol_gwas --subtype clin_gwas --Pcut 5e-8 --outfile "./output/cpag_output_molGWAS_vs_clinGWAS_pcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype clin_gwas --Pcut 5e-8 --outfile "./output/cpag_output_INTRA_clinGWAS_pcut5e-8_out.csv"

#python main.py cpagdb --threads 6 --subtype H2P2 --H2P2-Pcut 1e-5 --subtype NHGRI --outfile "./output/cpag_output_H2P2p1e5_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype NHGRI --outfile "./output/cpag_output_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 4 --subtype BloodMetabolites --subtype NHGRI --outfile "./output/cpag_output_BloodMetabolites_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 4 --subtype UrineMetabolites --subtype NHGRI --outfile "./output/cpag_output_UrineMetabolites_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 4 --subtype BloodXenobiotic --subtype NHGRI --outfile "./output/cpag_output_BloodXenobiotic_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype H2P2 --H2P2-Pcut 1e-4 --subtype NHGRI --outfile "./output/cpag_output_H2P2p1e4_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype H2P2 --H2P2-Pcut 1e-5 --subtype NHGRI --outfile "./output/cpag_output_H2P2p1e5_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype H2P2 --H2P2-Pcut 1e-6 --subtype NHGRI --outfile "./output/cpag_output_H2P2p1e6_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype H2P2 --H2P2-Pcut 1e-7 --subtype NHGRI --outfile "./output/cpag_output_H2P2p1e7_vs_NHGRIpcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype H2P2 --H2P2-Pcut 5e-8 --subtype NHGRI --outfile "./output/cpag_output_H2P2p5e8_vs_NHGRIpcut5e-8_out.csv"

#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_BloodMetabolites_vs_NHGRIpcut5e-8_out.csv
#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_UrineMetabolites_vs_NHGRIpcut5e-8_out.csv
#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_BloodXenobiotic_vs_NHGRIpcut5e-8_out.csv

#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_H2P2p1e5_vs_NHGRIpcut5e-8_out.csv
#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_H2P2p1e6_vs_NHGRIpcut5e-8_out.csv
#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_H2P2p1e7_vs_NHGRIpcut5e-8_out.csv
#python main.py post_analysis --anno-ontology --anno-cols 'Trait2' --infile ./output/cpag_output_H2P2p5e8_vs_NHGRIpcut5e-8_out.csv


## run covid19 #
# Ellingus data
#python main.py usr-gwas --threads 10 --infile infiles/covid19raw/covid19_ukbb/top_EllinghausPCs_covid19.csv --SNPcol "avsnp150" --delimitor "," --Pcol "p_value" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --ld-clump 1 --outfile output/top_EllinghausPCs_covid19_pcut1e-5.cpag2_out_20201022.csv



### run CPAG1
#python main.py cpagdb --threads 6 --subtype "CPAG1" --Pcut 5e-8 --outfile "./output/cpag1_20130904_output_INTRA_NHGRI_pcut5e-8_out.csv"
#python main.py cpagdb --threads 6 --subtype "CPAG1" --Pcut 1e-7 --outfile "./output/cpag1_20130904_output_INTRA_NHGRI_pcut1e-7_out.csv"
python main.py cpagdb --threads 6 --subtype "CPAG1" --Pcut 1e-6 --outfile "./output/cpag1_20130904_output_INTRA_NHGRI_pcut1e-6_out.csv"
python main.py cpagdb --threads 6 --subtype "CPAG1" --Pcut 1e-5 --outfile "./output/cpag1_20130904_output_INTRA_NHGRI_pcut1e-5_out.csv"



## run covid19
## 1) from covid19.hg
## 2) from NHBLI 

#for i in infiles/covid19raw/covid19_round4/*txt
for i in infiles/covid19raw/covid19_round4/COVID19_HGI_D1_ALL_20201020.b37_1.0E-5.txt
do
	echo $i
	#continue
	#python main.py usr-gwas --threads 3 --infile $i --SNPcol "rsid" --Pcol "all_inv_var_meta_p" --ld-clump 1 --usr-pcut 5e-8 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --outfile "$i"pcut5e8.cpag2_out.csv
	
done









## directly from manuscript supp. Table 5
#python main.py usr-gwas --threads 10 --infile "infiles/covid19raw/covid_19_table2.csv" --SNPcol "rsID" --delimitor "," --Pcol "p" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --outfile "./output/cpag_output_covid19table2old_suggestivepcut1e-5_out.csv"
#python main.py usr-gwas --threads 10 --infile "infiles/covid19raw/covid_19_table2_broad.csv" --ld-clump 1 --SNPcol "rsID" --delimitor "," --Pcol "pval" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --outfile "./output/cpag_output_covid19table2broad_suggestivepcut1e-5_out.csv"
#python main.py usr-gwas --threads 10 --infile "infiles/covid19raw/covid_19_table2_narrow.csv" --ld-clump 1  --SNPcol "rsID" --delimitor "," --Pcol "pval" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --outfile "./output/cpag_output_covid19table2narrow_suggestivepcut1e-5_out.csv"



## 2) from ukbb (NHGLI)
#python main.py usr-gwas --threads 10 --infile infiles/covid19raw/covid19_ukbb/top_EllinghausPCs_covid19.csv --SNPcol "avsnp150" --delimitor "," --Pcol "p_value" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --ld-clump 1 --outfile outfile/top_EllinghausPCs_covid19.txt_H2P2LDdb.cpag2_out_WithLDclump.csv

#python main.py usr-gwas --cpu 10 --infile infiles/covid19raw/covid19_ukbb/top_UKBB_covid19_ALLtested_061820.txt.gz --SNPcol "rsid" --delimitor "\t" --Pcol "p.value" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --H2P2-Pcut 1e-5 --ld-clump 1 --outfile infiles/covid19raw/covid19_ukbb/top_UKBB_covid19_ALLtested_061820.txt_1KGLDdb.cpag2_out.csv

#python main.py usr-gwas --cpu 10 --infile infiles/covid19raw/covid19_ukbb/top_UKBB_covid19_ALLtested_061820.txt.gz --SNPcol "rsid" --delimitor "\t" --Pcol "p.value" --usr-pcut 1e-5  --ld-clump 1 --outfile infiles/covid19raw/covid19_ukbb/top_UKBB_covid19_ALLtested_061820.txt_H2P2GLDdb.cpag2_out_NoLDclump.csv


#python main.py usr-gwas --cpu 10 --infile infiles/covid19raw/covid19_ukbb/top_UKBB_covid19_ALLtested_061820_lw.txt.gz --SNPcol "rsID_snpnexus" --delimitor "\t" --Pcol "p.value" --usr-pcut 1e-5 --cpagdb-pcut 5e-8 --ld-clump 1 --ld-clump-r2 0.6 --H2P2-Pcut 1e-5 --outfile infiles/covid19raw/covid19_ukbb/top_UKBB_covid19_ALLtested_061820.txt_snpnexusAnnot_1KGLDdb.cpag2_out_withLDclump.csv


	
	

#chmod g+r GWAS-QTL_data_dictionary.xlsx
scp -p GWAS-QTL_data_dictionary.xlsx minerva:/sc/arion/projects/ad-omics/data/references/GWAS/
ssh minerva "chmod g+wr /sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"


## METASOFT - meta-analysis of eQTLs

## Steps:

*createInputFile*

For a selection of datasets, extract the SNP, phenotype, beta and se and merge together into a single file with the following columns:

SNP-Phenotype pair, Beta1, SE1, Beta2, SE2, Beta3,SE3, .... BetaN, SE_N for N tissues

Read in by chromosome or by chunk to keep memory requirements reasonable

*runMETASOFT*

run METASOFT on the input file


# Script for Fine-mapping QTL and GWAS SNPs simulataneously
options(echo=TRUE)
# input - a COLOC results file
# RData object containing nested lists:
## locus_name
## .. gene_name
## .... COLOC summary
## ...... object
## ......... summary again
## ......... results ( coloc inputs for both GWAS and QTL)

# we're only interested in GWAS-QTL pairs that colocalise at PP4 > 0.5

library(tidyverse)
library(echolocatoR)

# for each locus, find the genes with PP4 > 0.5
# output - just the genes with a PP4 > threshold
prep_locus <- function(locus, pp4_threshold = 0.5){
    pp4 <- purrr::map_df( locus ,~{.x$df}, .id = "gene" ) %>% dplyr::filter(PP.H4.abf > pp4_threshold)
    if( nrow(pp4) == 0){ return(NULL) }
    locus <- locus[ pp4$gene ]
    return(locus)
}

# for each gene, prepare data for echoR

#inputfile <- "Microglia_all_regions_Kunkle_2019_COLOC.RData"
#outFolder <- "/hpc/users/humphj04/pipelines/downstream-QTL/Fine_Mapping/Kunkle_Microglia_all_regions"

# MyND monocytes - larger sample size = better fine-mapping?
inputfile <- "Monocytes_MyND_Kunkle_2019_COLOC.RData"
outFolder <- "/hpc/users/humphj04/pipelines/downstream-QTL/Fine_Mapping/Kunkle_Monocytes_MyND/"

# microglia young
#inputfile <- "Microglia_Young_Kunkle_2019_COLOC.RData"
#outFolder <- "/hpc/users/humphj04/pipelines/downstream-QTL/Fine_Mapping/Kunkle_Microglia_Young/"


if(!dir.exists(outFolder) ){
    dir.create(outFolder)
}

# name is all_obj
load(inputfile)

# extract colocalized genes
prepped <- map(all_obj, prep_locus)

# just use BIN1 for now
#prepped <- list(prepped$BIN1)
#prepped <- prepped[1:10]

# use Brian's function to write objects out as echoR input dataframes
merge_coloc_results(prepped, results_level = c("snp"), save_path = outFolder)

print("conversion to echoR completed")

#exit()
#stop()
# echoR 

# unzip the input files first!

fullSS_path <- file.path(outFolder, "merged_coloc_results.snp.tsv.gz")

top_SNPs <- import_topSNPs(
  topSS = fullSS_path,
  chrom_col = "chr",
  position_col = "pos",
  snp_col="snp",
  pval_col="gwas.pvalues",
  effect_col="gwas.beta",
  gene_col="gene",
  locus_col = "Locus",
  grouping_vars = c("Locus","Gene"))

loci <- gene_locus_list(top_SNPs)

# run the big function
gwas_res <- paste0(outFolder, "/GWAS_finemap_res.tsv")

if(FALSE){
#if( !file.exists(gwas_res) ){

finemap_res <- finemap_loci(# GENERAL ARGUMENTS 
                                        top_SNPs = top_SNPs,
                                        results_dir = outFolder,

                                        loci = loci,
                                        dataset_name = "Kunkle_GWAS",
                                        dataset_type = "GWAS",
                                        force_new_subset = TRUE,
                                        force_new_LD = F,
                                        force_new_finemap = TRUE,
                                        remove_tmps = FALSE,

                 
                 # SUMMARY STATS ARGUMENTS
                 fullSS_path = fullSS_path,
                 query_by ="coordinates",
                 chrom_col = "chr",
                 position_col = "pos",
                 snp_col = "snp",
                 pval_col = "gwas.pvalues",
                 effect_col = "gwas.beta",
                 stderr_col = "gwas.se",
                 freq_col = "gwas.MAF",
                 MAF_col = "calculate",
                 gene_col = "gene",
                 N_cases_col = "N_cases",
                 N_controls_col = "N_controls",
                 A1_col = "A1",
                 A2_col = "A2",
                 
                 # QTL prefixes
                 QTL_prefixes = c("qtl."),

                 # FILTERING ARGUMENTS
                 bp_distance = 500000, #100000,
                 min_MAF = 0.001,

                 # FINE-MAPPING ARGUMENTS
                 finemap_methods = c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                 n_causal = 5,
                 PP_threshold = .95,

                 # LD ARGUMENTS 
                 LD_reference = "1KGphase3",
                 superpopulation = "EUR",
                 download_method = "axel",

                 # PLOT ARGUMENTS 
                 plot.types=c("simple"),
                 plot.window = 500000/2,
                 
                 # CONDA
                 conda_env = "echolocatoR",   
                 )

# write out finemap res
readr::write_tsv(finemap_res, path = gwas_res )

}

# now for QTLs!

qtl_res <- paste0(outFolder, "/QTL_finemap_res.tsv")
if( TRUE){


fullSS_path <- file.path(outFolder, "merged_coloc_results.snp.tsv.gz")

top_SNPs <- import_topSNPs(
  topSS = fullSS_path,
  chrom_col = "chr",
  position_col = "pos",
  snp_col="snp",
  pval_col="qtl.pvalues",
  effect_col="qtl.beta",
  gene_col="gene",
  locus_col = "gene",
  grouping_vars = c("Locus","Gene"))

loci <- gene_locus_list(top_SNPs)


#if( !file.exists(qtl_res) ){
finemap_res <- finemap_loci(# GENERAL ARGUMENTS 
                                        top_SNPs = top_SNPs,
                                        results_dir = outFolder,

                                        loci = loci,
                                        dataset_name = "Kunkle_QTL",
                                        dataset_type = "QTL",
                                        force_new_subset = TRUE,
                                        force_new_LD = TRUE,
                                        force_new_finemap = TRUE,
                                        remove_tmps = FALSE,

                 # SUMMARY STATS ARGUMENTS
                 fullSS_path = fullSS_path,
                 query_by ="coordinates",
                 chrom_col = "chr",
                 position_col = "pos",
                 snp_col = "snp",
                 pval_col = "qtl.pvalues",
                 effect_col = "qtl.beta",
                 stderr_col = "qtl.se",
                 freq_col = "qtl.MAF",
                 MAF_col = "calculate",
                 sample_size = 90, # different for each QTL dataset
                 gene_col = "gene",
                 A1_col = "A1",
                 A2_col = "A2",

                 # QTL prefixes
                 QTL_prefixes = c("qtl."),

                 # FILTERING ARGUMENTS
                 bp_distance = 500000, #100000,
                 min_MAF = 0.001,

                 # FINE-MAPPING ARGUMENTS
                 finemap_methods = c("ABF","FINEMAP","SUSIE","POLYFUN_SUSIE"),
                 n_causal = 5,
                 PP_threshold = .95,

                 # LD ARGUMENTS 
                 LD_reference = "1KGphase3",
                 superpopulation = "EUR",
                 download_method = "axel",

                 # PLOT ARGUMENTS 
                 plot.types=c("simple"),
                 plot.window = 500000/2,

                 # CONDA
                 conda_env = "echolocatoR",
                 )

# write out finemap res
readr::write_tsv(finemap_res, path = qtl_res )
}


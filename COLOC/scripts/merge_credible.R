# input:

# directory containing GWAS files
# BED file of GWAS hits to test
# directory containing QTL nominal results 

# process:
# take the GWAS hit coordinates
# expand by 1MB either side
# extract the GWAS P-value distribution using tabix

# read in QTL parquet file for that chromosome
# extract SNP-Gene pairs at the coordinates
# split pairs by Gene

# for each gene:
# get together the two distributions, merge by shared position
# run COLOC to get different colocalisations
# extract top QTL SNP for locus
# record results in table - GWAS locus, Gene, top QTL SNP, COLOC probabilities
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(optparse))

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))


#library(arrow)

# flank coordinates by set number of bases (default 1MB)
# work on either a coordinate string or a dataframe containing chr start and end columns
option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = ""),
        make_option(c('--debug'), help = "load all files and then save RData without running COLOC", action = "store_true", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFolder <- opt$outFolder

debug <- opt$debug
#gwas_prefix <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/processed/Nicolas_2018_processed_"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_expression/peer30/LumbarSpinalCord_expression_peer30"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_splicing/peer20/LumbarSpinalCord_splicing_peer20"

#hits_file <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/COLOC/Nicolas_2018/Nicolas_2018_hits_1e-7.tsv"
#outFile <- "test_coloc_results.tsv"

#options(echo = TRUE)
# outFolder <- '/sc/arion/projects/bigbrain/data/ROSMAP/analysis/downstream-QTL/COLOC/'

main <- function(){
    
        #gwas_dataset <- "Kunkle_2019"
    #gwas_dataset <- "Nalls23andMe_2019"
    #gwas_dataset <- "Marioni_2018"
    #qtl_dataset <- "Microglia_THA"
    
    
    all_files <- list.files(paste0(outFolder), pattern = "_cre.tsv",
                            recursive = TRUE,
                            full.names = TRUE )
  
    print(all_files)
    
    read_cre <- function(filename){
      res <- read_tsv(filename) 
      if(nrow(res)==0){return(NULL)}
      res$A1 <- as.character(res$A1)
      return(res)
      
    }
      
    all_credible <- purrr::map_df(all_files, ~{
      read_cre(.x) 
      #dplyr::filter( nsnp > SNP_threshold ) #%>%
      #mutate(locus = as.character(locus)) %>%
      #mutate( QTL_chr = fix_chr(QTL_chr), GWAS_chr = fix_chr(GWAS_chr), GWAS_P = as.numeric(GWAS_P) )
    },)
    outFile <- paste0(outFolder, "all_COLOC_credible_results_merged_H4_0.5.tsv.gz")
    readr::write_tsv(all_credible, path = outFile)
}

main()


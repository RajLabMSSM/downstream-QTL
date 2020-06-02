# merge all COLOC results
# subset to all results with H4 >0.5
library(tidyverse)
library(LDlinkR)

inFolder <- "/sc/hydra/projects/ad-omics/microglia_omics/COLOC/"

all_files <- list.files(inFolder, pattern = "COLOC.tsv", recursive = TRUE, full.names = TRUE )

names(all_files) <- all_files

all_res <- purrr::map_df(all_files, ~{
    read_tsv(.x) %>% 
    filter( PP.H4.abf > 0.5 ) %>% 
    mutate(locus = as.character(locus)) 
}, .id = "file")

all_res <- mutate(all_res, file = gsub("_COLOC.tsv", "", basename(file) ) )


# get out QTL and GWAS names

# deal with Nicolas - only GWAS with 3 strings, rest are 2 strings
all_res$file <- gsub("_hg38", "", all_res$file)

all_res$QTL <- map_chr(str_split(all_res$file, "_"), ~{ n = length(.x); paste0( .x[1:(n-2)], collapse = "_") })

all_res$GWAS <- map_chr(str_split(all_res$file, "_"), ~{ n = length(.x); paste0(.x[n-1], "_", .x[n] ) })

# match in disease from GWAS name
gwas_key <- 
tribble(
  ~GWAS, ~disease,
    "Ripke_2014", "SCZ",
    "Wray_2018", "MDD",
    "IMSGC_2019", "MS",  
    "Stahl_2019", "BPD",
    "Daner_2020", "BPD",
    "Nalls23andMe_2019", "PD",
    "Lambert_2013", "AD",
    "Marioni_2018", "AD",
    "Jansen_2018",  "AD",
    "Kunkle_2019",  "AD",
    "Nicolas_2018", "ALS",
)

all_res <- left_join(all_res, gwas_key, by = "GWAS")

all_res <- select(all_res, GWAS, disease, QTL, -file, everything() )

# add pairwise LD using LDlink
# split into chunks by GWAS SNP
all_snps <- select(all_res, QTL_SNP, GWAS_SNP) %>% distinct() %>% split(.$GWAS_SNP)
# for each chunk generate LD matrix

calc_LD <- function( x ){
    print(x)
    snps <- unique( c(x$GWAS_SNP, x$QTL_SNP) )
    # get LDlink token from .Renviron (only Jack has this, for access go to https://ldlink.nci.nih.gov/?tab=apiaccess)
    token <- Sys.getenv("LDLINK_TOKEN")
    if( token == ""){
        warning(" * LDlink token not found!" )
        return(NA)
    }
    
    # if only one SNP then QTL and GWAS must be same SNP
    if( length(snps) == 1){ 
        x$LD <- 1; return(x) 
    }
    
    # get pairwise LD matrix
    ld_matrix <-  LDmatrix( snps = snps, pop = "CEU", r2d = "r2", token = Sys.getenv("LDLINK_TOKEN") )

    #stopifnot( nrow(ld_matrix) < length(snps) )

    ld_matrix <- tibble::column_to_rownames(ld_matrix, var = "RS_number")
    
    # if GWAS SNP isn't found in LDlink
    gwas_snp <- unique(x$GWAS_SNP)
    if( !gwas_snp %in% colnames(ld_matrix) ){
        x$LD <- NA
        return(x)
    }

    ld_select <-  ld_matrix[ x$QTL_SNP, x$GWAS_SNP]
    if( length(ld_select) > 1){
        ld_select <- ld_select[,1]
    }   

    # weird error when LDlink can't match all SNPs
    if( is.null(ld_select) ){
        x$LD <- NA; return(x)
    }
    stopifnot(nrow(x) == length(ld_select) )
    x$LD <- ld_select
    # wait 5 seconds before returning - makes sure API queries are spread out
    #Sys.sleep(time = 5)
    return(x)
}

all_snps_ld <- map_df(all_snps, calc_LD )

all_res <- left_join(all_res, all_snps_ld, by = c("GWAS_SNP", "QTL_SNP") )

# write out
write_tsv(all_res, paste0(inFolder, path = "all_COLOC_results_merged_H4_0.5.tsv") )

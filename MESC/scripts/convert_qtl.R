
# Convert mmQTL summary statistics to MESC format
# Jack Humphrey 2023

# Format:
#GENE: Gene name
#GENE_COORD: Chromosomal coordinate of gene midpoint.
#SNP: SNP ID
#CHR: Chromosome number of SNP
#SNP_COORD: Chromosomal coordinate of SNP
#N: Sample size of eQTL study (may differ by SNP)
#Z: eQTL Z-score
# 
# Files must be sorted by chromosome and then gene.
# GENE_COORD needs to be added in from metadata 
# N - this is fiddly, for each association we need to know how many cohorts were used. User must provide sizes for each cohort.

# inputs
# SUMSTATS file 
# COHORT SIZES - for N
# FEATURE METADATA - for coordinates

# read in sum stats - use tabix to chunk
# if QTL nominal stats stored in single tabixed file

library(optparse)
library(tidyverse)
option_list <- list(
    make_option(c('-q', '--QTL' ), help='The ID of the QTL dataset to use, from the database', default = ""),
    #make_option(c('-m', '--meta'), help = 'The path to the metadata file with coordinates for each feature', default = ""),
    #make_option(c('-s', '--sizes'), help = 'a TSV file of cohort sizes in same order as meta-analysis', default = "" ),
    make_option(c('-c', '--chr'), help = 'which chromosome to process', default = "21"),
    make_option(c('-p', '--prefix'), help = "output prefix") 
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

qtl_name <- opt$QTL
#meta_file <- opt$meta
#size_file <- opt$sizes
chr <- opt$chr
prefix <- opt$prefix

stopifnot(chr %in% 1:22)
out_file <- paste0(prefix, ".", chr, ".input.tsv")


pullData <- function(dataset, type = "GWAS"){
    message(Sys.time()," * selected dataset: ", dataset)
    db_path <- "/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"

    message(Sys.time()," * reading GWAS database from ", db_path)
    stopifnot( file.exists(db_path) )

    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = type, na= c("", "-","NA")))

    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )


    #stopifnot( file.exists(gwas$full_path) )

    return(gwas)
}


extractQTL_tabix <- function(qtl, coord){
    if( qtl$full_chrom_type == "chr1" ){
        if( !grepl("chr", coord) ){
            coord <- paste0("chr", coord)
        }
    }
    stopifnot( file.exists(qtl$full_path) )
    cmd <- paste( "ml bcftools; tabix -h ", qtl$full_path, coord )
    message("       * running command: ", cmd)
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    # deal with regions of no QTL association
    if( nrow(result) == 0){
        return(NULL)
    }
    #stopifnot( nrow(result) > 0 )
    message("       * QTL extracted!")
    return(result)
}


extractQTL <- function(qtl, coord, sig_level = 0.05, force_maf = FALSE){
    if( qtl$full_file_type == "tabix" ){
        result <- extractQTL_tabix(qtl, coord)
    }
    # if no SNPs found in region return NULL value
    if(is.null(result) ){ return(NULL) }
    #return(result)
    # assign column names
    cmd <- paste( "zless ", qtl$full_path, " | head -1 " )
    columns <- colnames(data.table::fread( cmd = cmd ))
    col_dict <- setNames(1:length(columns), columns)
    colnames(result) <- columns
    return(result)
}

# get count of samples used for each association
getCohortSize <- function(res, cohort_df){
    # use SEs from each cohort - if 0 assume cohort was not present
    df <- dplyr::select(res, contains("sd_tissue") )
    df <- df > 0
    x <- rowSums(sweep(df, MARGIN = 2, STATS = cohort_df$size, FUN = "*" ), na.rm = TRUE)
    res$N <- x
    return(res)
}

matchCoords <- function(res, meta_df){
    meta_df$coord <- meta_df$start + floor((meta_df$end - meta_df$start)/2)
    
    # deal with sQTLs
    if( sum(grepl(":", res$feature)) > 0 ){
        if( sum(grepl("clu_", res$feature) > 0) ){
            res$gene <- str_split_fixed(res$feature, ":", 5)[,5]
        }else{
            res$gene <- str_split_fixed(res$feature, ":", 4)[,4]
        }
        res$GENE_COORD <- meta_df$coord[match(res$gene, meta_df$feature)]
    }else{

        # edge case - some sumstats have ensembl IDs without tag numbers
        if( sum(res$feature %in% meta_df$feature ) == 0 ){
            meta_df$feature <- gsub("\\.[0-9]+","", meta_df$feature )
        }
        res$GENE_COORD <- meta_df$coord[match(res$feature, meta_df$feature)]
    }
    return(res)
}

standardiseCols <- function(res, qtl){
    names(res)[ names(res) == qtl$full_pheno ] <- "feature"
    names(res)[ names(res) == qtl$full_snp ] <- "variant_id"   
    names(res)[ names(res) == qtl$full_chrom ] <- "chr"
    names(res)[ names(res) == qtl$full_pos ] <- "pos"
    names(res)[ names(res) == qtl$full_effect ] <- "beta"
    names(res)[ names(res) == qtl$full_se ] <- "se"
    # weird edge case where there are two separate chr columns with different names
    if( sum(names(res) == "chr" ) > 1 ){
        res <- res[, -which(names(res) == "chr")[1] ]
    }
    return(res)
}

calcZ <- function(res){
    if("Random_Z" %in% names(res) ){
        res$Z <- res$Random_Z
    }else{
        res$Z <- res$beta / res$se 
    }
    return(res)
}

#test_coord <- "chr21"
#qtl_name <- "mmQTL_GENCODE_expression"
#meta_file <- "/hpc/users/humphj04/pipelines/downstream-QTL/MESC/example/gencode.v38.gene_pheno_meta.tsv"
#cohort_file <- "/hpc/users/humphj04/pipelines/downstream-QTL/MESC/example/test_cohort_sizes.txt"

#out_file <- "test_mesc_input.tsv"

qtl <- pullData(qtl_name, type = "QTL")
# load in QTLs
res <- extractQTL(qtl, coord = chr)
print(head(res))
# standardise column names 
res <- standardiseCols(res,qtl)
print(head(res))

## feature coordinates - use GENCODE v30 gene coordinates by default
if( is.na(qtl$full_feature_meta) ){
    meta_file <- "/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.genes.tsv"
}else{
    meta_file <- qtl$full_feature_meta
}
stopifnot(file.exists(meta_file))

# if QTL comes from meta-analysis then you require a cohort sizes file
# otherwise use N provided
if( !is.na(qtl$full_cohort_meta) ){

    size_file <- qtl$full_cohort_meta
    stopifnot(file.exists(size_file))
    cohort_df <- readr::read_tsv(size_file)

    res <- getCohortSize(res,cohort_df)
}else{
    stopifnot(!is.na(qtl$N ) )
    res$N <- qtl$N
}

print(head(res))

meta_df <- readr::read_tsv(meta_file)

res <- res  %>%
    matchCoords(meta_df) 
print(head(res))
res <- res %>%
    calcZ() 
print(head(res))

res <- res %>%
    select(GENE = feature, GENE_COORD, SNP = variant_id, CHR = chr, SNP_COORD = pos, N, Z ) %>%
    arrange(GENE, SNP_COORD) %>%
    mutate(CHR = gsub("chr", "", CHR) ) %>%
    drop_na()

stopifnot( nrow(res) > 0)
print(head(res))


write_tsv(res, out_file)


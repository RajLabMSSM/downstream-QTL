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
#library(echolocatoR)
library(optparse)


#inputfile <- "Microglia_all_regions_Kunkle_2019_COLOC.RData"
#outFolder <- "Microglia_all_regions_Kunkle_2019/"
#QTL <- "Microglia_all_regions"
#GWAS <- "Kunkle_2019"


option_list <- list(
    make_option(c('--inFile'), help='', default = "example"),
    make_option(c('--outFolder'), help='', default = "example"),
    make_option(c('-g', '--GWAS'), help='', default = "Kunkle_2019"),
    make_option(c('-q','--QTL'), help='', default = "Microglia_all_regions") 
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

inFile <- opt$inFile
outFolder <- opt$outFolder
GWAS <- opt$GWAS
QTL <- opt$QTL


dataCode <- paste0(QTL, "_", GWAS)
#outFolder <- file.path(outFolder, dataCode)

outFile <-  paste0(outFolder, "/", dataCode,"_coloc_results.snp.1kgp3_ld.tsv.gz") 


# GWAS 
gwas_folder <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/"
stopifnot(dir.exists(gwas_folder) )

available_gwas <- list.files(path = gwas_folder, include.dirs=TRUE)

stopifnot(GWAS %in% available_gwas)


# replace gene IDs with gene names
gene_meta <- readr::read_tsv("/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.tx2gene.tsv") %>%
    select(GENENAME, GENEID) %>% distinct() %>%
    mutate(GENEID_notag = gsub(".[0-9]+$", "", GENEID) )

# for each locus, find the genes with PP4 > 0.5
# output - just the genes with a PP4 > threshold
prep_locus <- function(locus, pp4_threshold = 0.5){
    pp4 <- purrr::map_df( locus ,~{.x$df}, .id = "gene" ) %>%
        dplyr::filter(PP.H4.abf > pp4_threshold)
    if( nrow(pp4) == 0){ return(NULL) }
    locus <- locus[ pp4$gene ]
    return(locus)
}

map_gene <- function(locus, tag = FALSE){
    if( tag == TRUE){
        geneid_vector <- gene_meta$GENEID
    }else{
        gene_id_vector <- gene_meta$GENEID_notag
    }
    # for each QTL gene in locus replace gene id with gene name
    # replace locus with paste0(locus, ":", genename)
    if( is.null(locus) ){ return(NULL) }
    names(locus) <- gene_meta$GENENAME[ match(names(locus), gene_id_vector) ]
    locus <- map(locus, ~{
        gene_df <- .x$object$results
        gene_df$gene <- gene_meta$GENENAME[ match(gene_df$gene, gene_id_vector) ]
        .x$object$results <- gene_df
        return(.x)
    })
    return(locus)
}




# for each gene, prepare data for echoR

if(!dir.exists(outFolder) ){
    dir.create(outFolder, recursive = TRUE)
}

# name is all_obj
load(inFile)

# extract colocalized genes
prepped <- map(all_obj, prep_locus) %>% purrr::compact()

# convert gene id to gene name
# drop for now - breaks on sQTLs
#prepped <- map(prepped, map_gene)

save(prepped, file =  "test.RData")

# use Brian's function to write objects out as echoR input dataframes
echolocatoR::merge_coloc_results(prepped, results_level = c("snp"), save_path = outFolder)

print("conversion to echoR completed")

### Get LD for each GWAS locus

merged_loci <- read_tsv( paste0(outFolder, "/merged_coloc_results.snp.tsv.gz"))

## Harmonise locus names with Brian's fine-mapping
refFolder <- "/sc/arion/projects/ad-omics/data/references/GWAS/"

## AD - Rename Loci to be standard
if( GWAS %in% c("Lambert_2013", "Marioni_2018", "Jansen_2018", "Kunkle_2019") ){
ad_consensus_loci <- readxl::read_excel(paste0(refFolder,"AD_GWAS_consensus_loci_names.xlsx") )

merged_loci <- left_join(merged_loci, ad_consensus_loci, by = c("Locus" = "consensus_locus_name")) 
merged_loci$Locus <- coalesce(merged_loci$old_locus_name, merged_loci$Locus)
}

## MS - Brian used different locus names in MS for fine-mapping, use a key
if( GWAS == "IMSGC_2019"){
ms_locus_key <- read_tsv(paste0(refFolder,"MS_GWAS_locus_key.txt"))

merged_loci$Locus_Brian <- ms_locus_key$Locus_Brian[ match(merged_loci$Locus, ms_locus_key$Locus)]
merged_loci$Locus <- coalesce(merged_loci$Locus_Brian, merged_loci$Locus)
merged_loci$Locus_Brian <- NULL
}

if( GWAS == "Stahl_2019"){
# Stahl - remove "*" and gsub "," to "-"
merged_loci$Locus <- gsub("\\*", "", gsub(",", "-", merged_loci$Locus))
}

if( GWAS == "Nalls23andMe_2019"){
# Nalls PD - 5 loci were skipped by Brian because they didn't pass QC
merged_loci <- filter(merged_loci, !Locus %in% c("DDRGK1", "SEMA4A", "GXYLT1", "FGD4", "ZNF608"))
}

# split merged loci back into individual loci
all_loci_split <- split(merged_loci, merged_loci$Locus)

all_loci_ld <- map2_df(.y = all_loci_split, .x= names(all_loci_split), ~{
    print(paste0(' * adding LD for Locus: ', .x))
    locus_path <- file.path(gwas_folder, GWAS, .x)
    #stopifnot(dir.exists(locus_path) )
    if( !dir.exists(locus_path) ){ return(NULL) }
    
    ld_file <- file.path(gwas_folder, GWAS,.x, "LD", paste0( .x , ".1KGphase3_LD.RDS")  )
    stopifnot(file.exists(ld_file))
    
    print(ld_file)
    ld_matrix <- readRDS(ld_file)
    
    lead_gwas_snp <- arrange(.y,  gwas.pvalues, desc(abs(gwas.beta)) ) %>% head(1) %>% pull(snp)
    lead_qtl_snp <-  arrange(.y,  qtl.pvalues, desc(abs(qtl.beta)) ) %>% head(1) %>% pull(snp)
    
    print(paste0(" * lead QTL SNP: ", lead_qtl_snp) )
    print(paste0(" * lead GWAS SNP: ", lead_gwas_snp) )
    
    ld_res <- .y
    
    # stupid error handling - if either SNP is not in the LD matrix
    if(lead_qtl_snp %in% colnames(ld_matrix) ){
        lead_qtl_ld <- enframe(ld_matrix[,lead_qtl_snp])
        names(lead_qtl_ld) <- c("snp", "lead_qtl_ld")
        lead_qtl_ld$snp <- row.names(ld_matrix)
        
        ld_res <- left_join(.y, lead_qtl_ld, by = "snp")
        ld_res$lead_qtl_ld <- ld_res$lead_qtl_ld^2

    }else{
        ld_res$lead_qtl_snp <- NA
    }
    
    if(lead_gwas_snp %in% colnames(ld_matrix) ){
        lead_gwas_ld <- enframe(ld_matrix[,lead_gwas_snp])
        names(lead_gwas_ld) <- c("snp", "lead_gwas_ld")
        lead_gwas_ld$snp <- row.names(ld_matrix)    
        
        ld_res <- left_join(.y, lead_gwas_ld, by = "snp")
        ld_res$lead_gwas_ld <- ld_res$lead_gwas_ld^2

    }else{
        ld_res$lead_gwas_snp <- NA
    }
   
    #ld_res <- left_join(.y, lead_gwas_ld, by = "snp") %>% left_join(lead_qtl_ld, by = "snp")
    
    # the LD matrices store R, not R^2. So square them!
        
    return(ld_res)

})

write_tsv(all_loci_ld, path = outFile)

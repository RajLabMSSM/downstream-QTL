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

#inputfile <- "Microglia_all_regions_Kunkle_2019_COLOC.RData"
#outFolder <- "/hpc/users/humphj04/pipelines/downstream-QTL/Fine_Mapping/Kunkle_Microglia_all_regions"

# MyND monocytes - larger sample size = better fine-mapping?
inputfile <- "Microglia_all_regions_Nalls23andMe_2019_COLOC.RData"
#inputfile <- "Monocytes_MyND_Kunkle_2019_COLOC.RData"
outFolder <- "/hpc/users/humphj04/pipelines/downstream-QTL/Fine_Mapping/Microglia_all_regions_Nalls23andMe_2019/"

inputfile <- "Microglia_all_regions_Kunkle_2019_COLOC.RData"
outFolder <- "Microglia_all_regions_Kunkle_2019/"
GWAS <- "Kunkle_2019"

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

# convert gene id to gene name
prepped <- map(prepped, map_gene)

# just use BIN1 for now
#prepped <- list(prepped$BIN1)
#prepped <- prepped[1:10]

# use Brian's function to write objects out as echoR input dataframes
merge_coloc_results(prepped, results_level = c("snp"), save_path = outFolder)

print("conversion to echoR completed")

### Get LD for each GWAS locus

merged_loci <- read_tsv( paste0(outFolder, "/merged_coloc_results.snp.tsv.gz"))

## Harmonise locus names with Brian's fine-mapping

## Rename AD Loci
refFolder <- "/sc/arion/projects/ad-omics/data/references/GWAS/"

ad_consensus_loci <- readxl::read_excel(paste0(refFolder,"AD_GWAS_consensus_loci_names.xlsx") )

merged_loci <- left_join(merged_loci, ad_consensus_loci, by = c("Locus" = "consensus_locus_name")) 
all_res$Locus <- coalesce(all_res$old_locus_name, all_res$Locus)

## MS - Brian used different locus names in MS for fine-mapping, use a key
ms_locus_key <- read_tsv(paste0(refFolder,"MS_GWAS_locus_key.txt"))

merged_loci$Locus_Brian <- ms_locus_key$Locus_Brian[ match(merged_loci$Locus, ms_locus_key$Locus)]
merged_loci$Locus <- coalesce(merged_loci$Locus_Brian, merged_loci$Locus)
merged_loci$Locus_Brian <- NULL

# Stahl - remove "*" and gsub "," to "-"
merged_loci$Locus <- gsub("\\*", "", gsub(",", "-", merged_loci$Locus))

# Nalls PD - 5 loci were skipped by Brian because they didn't pass QC
merged_loci <- filter(merged_loci, !Locus %in% c("DDRGK1", "SEMA4A", "GXYLT1", "FGD4", "ZNF608"))


gwas_folder <- "/sc/arion/projects/pd-omics/brian/Fine_Mapping/Data/GWAS/"
stopifnot(dir.exists(gwas_folder) )

available_gwas <- list.files(path = gwas_folder, include.dirs=TRUE)

stopifnot(GWAS %in% available_gwas)

# split merged loci back into individual loci
all_loci_split <- split(merged_loci, merged_loci$Locus)

# for each locus find the lead QTL and GWAS SNPs
# extract LD for each SNP with both types of SNP
#.y <- all_loci_split$MED12L
#.x <- "MED12L"

all_loci_ld <- map2_df(.y = all_loci_split, .x= names(all_loci_split), ~{
    print(paste0(' * adding LD for Locus: ', .x))
    locus_path <- file.path(gwas_folder, GWAS, .x)
    #stopifnot(dir.exists(locus_path) )
    if( !dir.exists(locus_path) ){ return(NULL) }
    
    ld_file <- file.path(gwas_folder, GWAS,.x, "LD", paste0( .x , ".1KGphase3_LD.RDS")  )
    stopifnot(file.exists(ld_file))

    ld_matrix <- readRDS(ld_file)

    lead_gwas_snp <- arrange(.y,  gwas.pvalues, desc(abs(gwas.beta)) ) %>% head(1) %>% pull(snp)
    lead_qtl_snp <-  arrange(.y,  qtl.pvalues, desc(abs(qtl.beta)) ) %>% head(1) %>% pull(snp)

    lead_qtl_ld <- enframe(ld_matrix[,lead_qtl_snp])
    names(lead_qtl_ld) <- c("snp", "lead_qtl_ld")

    lead_gwas_ld <- enframe(ld_matrix[,lead_gwas_snp])
    names(lead_gwas_ld) <- c("snp", "lead_gwas_ld")

    ld_res <- left_join(.y, lead_gwas_ld, by = "snp") %>% left_join(lead_qtl_ld, by = "snp")
    
    # the LD matrices store R, not R^2. So square them!
    ld_res$lead_qtl_ld <- ld_res$lead_qtl_ld^2
    ld_res$lead_gwas_ld <- ld_res$lead_gwas_ld^2

    return(ld_res)

})

write_tsv(all_loci_ld, path = paste0(outFolder, GWAS, "_coloc_results.snp.1kgp3_ld.tsv.gz") )


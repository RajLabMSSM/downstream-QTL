# merge all COLOC results
# subset to all results with H4 >0.5
library(tidyverse)
library(LDlinkR)
library(optparse)
option_list <- list(
    make_option(c('--inFolder' ), help='The full path to the folder that contains the COLOC results', default = ""),
    make_option(c('--threshold'), help = "minimum PP.H4.abf value", default = 0),
    make_option(c('-g', '--geneMeta'), help = "Path to gene metadata, matching Ensembl IDs to gene names", default = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.tx2gene.tsv"),
    make_option(c('--ld'), help = "whether to match LD or not", action="store_true", default=FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

calculate_LD <- opt$ld
inFolder <- opt$inFolder
H4_threshold <- opt$threshold
geneMeta <- opt$geneMeta

if( calculate_LD == TRUE){
    LD_string <- "_with_LD"
}else{
    LD_string <- "_no_LD"
}

fix_chr <- function(chr_string){
    if( all(!grepl("chr", chr_string) ) ){
        chr_string <- paste0("chr", chr_string)
    }
    return(chr_string)
}

outFile <- paste0(inFolder, "all_COLOC_results_merged_H4_", H4_threshold, LD_string,".tsv.gz")

message(" * writing to ", outFile)

#inFolder <- "/sc/arion/projects/ad-omics/microglia_omics/COLOC/"
#H4_threshold <- 0.5
#calculate_LD <- TRUE
#H4_threshold <- 0
#calculate_LD <- FALSE

all_files <- list.files(inFolder, pattern = "COLOC.tsv", recursive = TRUE, full.names = TRUE )

names(all_files) <- all_files

all_res <- purrr::map_df(all_files, ~{
    read_tsv(.x) %>% 
    filter( PP.H4.abf >= H4_threshold ) %>% 
    mutate(locus = as.character(locus)) %>%
    mutate( QTL_chr = fix_chr(QTL_chr), GWAS_chr = fix_chr(GWAS_chr), GWAS_P = as.numeric(GWAS_P) )
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
    "Bellenguez_2021", "AD", 
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
    "NicolasSuggestive_2018", "ALS"
)

all_res <- left_join(all_res, gwas_key, by = "GWAS")

all_res <- select(all_res, GWAS, disease, QTL, -file, everything() )

# deal with gene
# sQTLs include junction too
all_res$feature <- all_res$gene

all_res$geneid <- map_chr(str_split(all_res$gene, ":"), ~{ .x[ length(.x) ] })


# TODO: if fusion or antisense gene then don't trim off tag
#all_res$geneid <- ifelse( grepl("_|AS", all_res$geneid), all_res$geneid,  str_split_fixed(all_res$geneid, "\\.", 2)[,1] )

# have separate junction column
all_res$QTL_junction <- map_chr(str_split(all_res$gene, ":"), ~{ paste0(.x[1], ":", .x[2], "-", .x[3])  })


gene_meta <- 
    read_tsv(geneMeta) %>%
#read_tsv("/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.tx2gene.tsv") %>% 
    janitor::clean_names() %>%
    select(genename, geneid) %>% distinct()

# remove tags
#gene_meta$geneid <- str_split_fixed(gene_meta$geneid, "\\.", 2)[,1]

# match on gene symbols
all_res$genename <- gene_meta$genename[ match(all_res$geneid, gene_meta$geneid) ]

# coalesce - if no gene symbol found use ID
all_res$QTL_Gene <- coalesce(all_res$genename, all_res$geneid)

# add geneid back to make sure
all_res$QTL_Ensembl <- gene_meta$geneid[match(all_res$QTL_Gene, gene_meta$genename)] 

#all_res$type <- ifelse( grepl("sQTL", all_res$QTL), "sQTL", "eQTL" )

# make junction NA if eQTL
#all_res$QTL_junction <- ifelse(all_res$type == "sQTL", all_res$QTL_junction, ".")

all_res <- select(all_res, disease, GWAS, locus, starts_with("GWAS"), QTL, feature, starts_with("QTL"), nsnps, starts_with("PP") )

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
   
    # weird snps with "esv" instead of "rs" - remove
    snps <- snps[ grepl("rs", snps) ]
    # also ignores non-rs IDs - in format 1:1:C:G
    # if only one SNP then QTL and GWAS must be same SNP
    if( length(snps) == 1){ 
        x$LD <- 1; return(x) 
    }
    
    # get pairwise LD matrix
    ld_matrix <-  LDmatrix( snps = snps, pop = "EUR", r2d = "r2", token = Sys.getenv("LDLINK_TOKEN") )

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
    Sys.sleep(time = 5)
    return(x)
}

if( calculate_LD == TRUE){
    all_snps_ld <- map_df(all_snps, calc_LD )

    all_res <- left_join(all_res, all_snps_ld, by = c("GWAS_SNP", "QTL_SNP") )
}

# write out
write_tsv(all_res, outFile)

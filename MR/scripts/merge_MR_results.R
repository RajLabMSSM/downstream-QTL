# merge all COLOC results
# subset to all results with H4 >0.5
library(tidyverse)
#library(LDlinkR)
library(optparse)
option_list <- list(
    make_option(c('--inFolder' ), help='The full path to the folder that contains the COLOC results', default = "./"),
    make_option(c('-t','--threshold'), help = "minimum SNP count ", default = 0),
    make_option(c('-f', '--fdr'), help = "the exposure filter option, if TRUE, use FDR < 0.05", default = FALSE),
    make_option(c('-g', '--geneMeta'), help = "Path to gene metadata, matching Ensembl IDs to gene names", default = "/sc/arion/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.tx2gene.tsv.gz")
    #make_option(c('--ld'), help = "whether to match LD or not", action="store_true", default=FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

#calculate_LD <- opt$ld
inFolder <- opt$inFolder
SNP_threshold <- opt$threshold
fdr_filter <- opt$fdr
geneMeta <- opt$geneMeta


# if( calculate_LD == TRUE){
#     LD_string <- "_with_LD"
# }else{
#     LD_string <- "_no_LD"
# }

fix_chr <- function(chr_string){
    if( all(!grepl("chr", chr_string) ) ){
        chr_string <- paste0("chr", chr_string)
    }
    return(chr_string)
}
#inFolder <- "/sc/arion/projects/bigbrain/data/ROSMAP/analysis/downstream-QTL/MR/"
#SNP_threshold <- 0
# fdr_filter <- FALSE
#fdr_filter <- TRUE
#SNP_threshold <- 1
# inFolder <- "/sc/arion/projects/bigbrain/data/ROSMAP/analysis/downstream-QTL/MR/"
# SNP_threshold <- 0
# fdr_filter <- TRUE
# fdr_filter <- FALSE

if(fdr_filter){
    outFile <- paste0(inFolder, "all_MR_results_merged_H4_", SNP_threshold, "_snp_wFDR.tsv.gz")
    message(" * writing to ", outFile)
    all_files <- list.files(inFolder, pattern = "_MR_wFDR.tsv", recursive = TRUE, full.names = TRUE )
}else{
    outFile <- paste0(inFolder, "all_MR_results_merged_H4_", SNP_threshold, "_snp.tsv.gz")
    message(" * writing to ", outFile)
    all_files <- list.files(inFolder, pattern = "_MR.tsv", recursive = TRUE, full.names = TRUE )
}


names(all_files) <- all_files

read_mr <- function(filename){
    res <- read_tsv(filename) 
    if(nrow(res)==0){return(NULL)}
    res$GWAS_P <- as.numeric(res$GWAS_P)
    return(res)
}

all_res <- purrr::map_df(all_files, ~{
    read_mr(.x)
    #dplyr::filter( nsnp > SNP_threshold ) #%>% 
    #mutate(locus = as.character(locus)) %>%
    #mutate( QTL_chr = fix_chr(QTL_chr), GWAS_chr = fix_chr(GWAS_chr), GWAS_P = as.numeric(GWAS_P) )
}, .id = "file") 
all_res["nsnp"][is.na(all_res["nsnp"])] <- 0

all_res <- all_res %>% dplyr::filter( nsnp > SNP_threshold ) %>% 
    dplyr::select(-id.outcome,-id.exposure,)
if(fdr_filter){
  all_res <- mutate(all_res, file = gsub("_MR_wFDR.tsv", "", basename(file) ) )
}else{
  all_res <- mutate(all_res, file = gsub("_MR.tsv", "", basename(file) ) )
}
all_res$exposure <- str_split_fixed(all_res$exposure, "\\.", 2)[,1]

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
"NicolasSuggestive_2018", "ALS",
"VanRheenenEUR_2021", "ALS",
"TrubetskoyEUR_2022", "SCZ",
"Yengo_2022", "Height",
"Ishigaki_2022", "RA"
)

all_res <- left_join(all_res, gwas_key, by = "GWAS")

all_res <- dplyr::select(all_res, GWAS, disease, QTL, -file, everything() )

# deal with gene
# sQTLs include junction too
all_res$feature <- all_res$gene
# 
# all_res$geneid <- map_chr(str_split(all_res$gene, ":"), ~{ .x[ length(.x) ] })


# TODO: if fusion or antisense gene then don't trim off tag
#all_res$geneid <- ifelse( grepl("_|AS", all_res$geneid), all_res$geneid,  str_split_fixed(all_res$geneid, "\\.", 2)[,1] )

# have separate junction column
#all_res$QTL_junction <- map_chr(str_split(all_res$gene, ":"), ~{ paste0(.x[1], ":", .x[2], "-", .x[3])  })


gene_meta <- 
    # read_tsv(geneMeta) %>%
read_tsv("/sc/arion/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.tx2gene.tsv.gz") %>%
    janitor::clean_names() %>%
    dplyr::select(genename, geneid) %>% distinct()

# remove tags
gene_meta$geneid <- str_split_fixed(gene_meta$geneid, "\\.", 2)[,1]

# match on gene symbols
all_res$exposure_Gene <- gene_meta$genename[ match(all_res$exposure, gene_meta$geneid) ]

# coalesce - if no gene symbol found use ID
all_res$exposure_Gene <- coalesce(all_res$exposure_Gene, all_res$exposure)

# add geneid back to make sure
all_res$exposure_Ensembl <- gene_meta$geneid[match(all_res$exposure_Gene, gene_meta$genename)] 

#colnames(all_res)
#all_res$type <- ifelse( grepl("sQTL", all_res$QTL), "sQTL", "eQTL" )

# make junction NA if eQTL
#all_res$QTL_junction <- ifelse(all_res$type == "sQTL", all_res$QTL_junction, ".")

all_res <- all_res%>% dplyr::select(-exposure)
all_res <- dplyr::select(all_res, disease, GWAS, locus, starts_with("GWAS"), QTL, feature, starts_with("QTL"), 
                         nsnp, outcome, starts_with("exposure"), method, nsnp, b, se, pval)
colnames(all_res)

head(all_res)
# add pairwise LD using LDlink
# split into chunks by GWAS SNP
# all_snps <- select(all_res, QTL_SNP, GWAS_SNP) %>% distinct() %>% split(.$GWAS_SNP)



# write out
message(" * writing to ", outFile)
write_tsv(all_res, outFile)

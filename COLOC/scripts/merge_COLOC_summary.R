# merge all COLOC results
# subset to all results with H4 >0.5
library(tidyverse)
library(LDlinkR)
library(optparse)
option_list <- list(
    make_option(c('-i', '--inFolder' ), help='The full path to the folder that contains the COLOC results', default = ""),
    make_option(c('-g', '--geneMeta'), help = "Path to gene metadata, matching Ensembl IDs to gene names", default = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.tx2gene.tsv")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

inFolder <- opt$inFolder
geneMeta <- opt$geneMeta
outFile <- paste0(inFolder, "all_COLOC_summary_results.tsv.gz")
message(" * writing to ", outFile)

all_files <- list.files(inFolder, pattern = "_COLOC_summary_level.tsv", recursive = TRUE, full.names = TRUE )

names(all_files) <- all_files

all_res <- purrr::map_df(all_files, ~{
    read_tsv(.x) %>% 
    mutate(locus = as.character(locus))
})


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

all_res <- select(all_res, GWAS, disease, QTL,  everything() )

# deal with gene
# sQTLs include junction too
#all_res$geneid <- map_chr(str_split(all_res$geneid, ":"), ~{ .x[ length(.x) ] })
all_res$geneid <- str_split_fixed(all_res$geneid, "\\.", 2)[,1]

# expects columns genename and geneid
#gene_meta <- read_tsv("/sc/hydra/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v30.tx2gene.tsv") %>% 
gene_meta <- read_tsv(geneMeta) %>%
    janitor::clean_names() %>%
    select(genename, geneid) %>% distinct()

# remove tags
gene_meta$geneid <- str_split_fixed(gene_meta$geneid, "\\.", 2)[,1]

# match on gene symbols
all_res$genename <- gene_meta$genename[ match(all_res$geneid, gene_meta$geneid) ]

# coalesce - if no gene symbol found use ID
all_res$QTL_Gene <- coalesce(all_res$genename, all_res$geneid)

# add geneid back to make sure
all_res$QTL_Ensembl <- gene_meta$geneid[match(all_res$QTL_Gene, gene_meta$genename)] 

all_res$type <- ifelse( grepl("sQTL", all_res$QTL), "sQTL", "eQTL" )

all_res <- select(all_res, disease, GWAS, locus, starts_with("GWAS"), QTL, type, starts_with("QTL"), nsnps, starts_with("PP") )

# write out
write_tsv(all_res, outFile)

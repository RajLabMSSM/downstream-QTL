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
gwas_key <- tribble(
  ~GWAS, ~disease,
  "Ripke_2014",              "SCZ",
  "Wray_2018",               "MDD",
  "IMSGC_2019",              "MS",
  "Stahl_2019",              "BPD",
  "Daner_2020",              "BPD",
  "Nalls23andMe_2019",       "PD",
  "Lambert_2013",            "AD",
  "Marioni_2018",            "AD",
  "Jansen_2018",             "AD",
  "Kunkle_2019",             "AD",
  "NicolasSuggestive_2018",  "ALS",
  "Pottier_2019",            "FTD",
  "Bellenguez_2021",         "AD",
  "VanRheenenEUR_2021",      "ALS",
  "VanRheenenASN_2021",      "ALS",
  "VanRheenenALL_2021",      "ALS",
  "Mullins_2021",            "BPD",
  "Farrell_2024",            "PSP",
  "TrubetskoyEUR_2022",      "SCZ",
  "LBD_AMPPD2023",           "LBD",
  "ABETA_2024",              "Abeta",
  "Yengo_2022",              "Height",
  "Ishigaki_2022",           "RA",
  "Belloy_2023APOE4",        "AD_APOE4",
  "Belloy_2023nonAPOE4",     "AD_noAPOE4",
  "Jones_2020",              "sCJD",
  "Chia_2024additive",       "MSA",
  "Chia_2024recessive",      "MSA",
  "Manzoni_2024",            "sFTLD",
  "MDD_2025",                "MDD",
  "Hatzikotoulas_2025",      "OA",
  "GP2_2025",                "PD",
  "GP2_clinical2025",        "PD_case_control",
  "Belloy_2025APOE4",        "AD_APOE4",
  "Belloy_2025nonAPOE4",     "AD_noAPOE4",
  "Chang_2017",              "PD"
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

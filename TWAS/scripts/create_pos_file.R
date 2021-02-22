# Create POS file

# this is for public data where they didn't provide a POS file
# 
library(tidyverse)
library(optparse)
library(GenomicRanges)
# inFolder
option_list <- list(
    make_option(c('-p', '--panel' ), help="the name of the panel - the data_code", default = ""),
    make_option(c('-w', '--weight_folder' ), help = "the path to the folder where the weight RDat files are kept", default = ""),
    make_option(c('-n', '--n_samples'), help = "the sample size for the QTLs")
)   

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

panel <- opt$panel
n_samples <- opt$n_samples
weight_folder <- opt$weight_folder

outFile <- paste0(weight_folder, "/", panel, ".pos")

weights <- list.files(weight_folder, pattern = ".RDat", full.names = TRUE)
#weights <- gsub("//", "/", weights)
#weights <- gsub(weight_folder, "", weights)

# assume ensembl id for now
genes <- gsub(paste0(panel, "."), "", gsub(".wgt.RDat", "", basename(weights) ) )

# read in GENCODE key for ensembl to gene name
#gtf <- read_tsv("/sc/arion/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v30.tx2gene.tsv") %>%
#    select(-TXNAME) %>% distinct()
gtf <- rtracklayer::import("/sc/arion/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v32.genes.gtf")

gtf$gene_id_no_tag <- gsub("\\.[0-9]+", "", gtf$gene_id)

if( sum(genes %in% gtf$gene_id_no_tag) > sum(genes %in% gtf$gene_id) ){
    gene_names <- gtf$gene_name[match(genes, gtf$gene_id_no_tag)]
}else{
    gene_names <- gtf$gene_name[match(genes, gtf$gene_id)]
}

df <- tibble(
  PANEL = panel,
  WGT = weights,
  ID = gene_names,
  )
df <- df[!is.na(df$ID),]

gtf_df <- gtf[ match(df$ID, gtf$gene_name) ]

df$CHR <- as.numeric(seqnames(gtf_df))
df$P0 <- start(gtf_df) - 1e6
df$P1 <- end(gtf_df) + 1e6
df$P0 <- ifelse(df$P0 < 0, 1, df$P0)
df$N <- n_samples

write_tsv(df, outFile)

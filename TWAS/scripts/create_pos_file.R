# Create POS file

# this is for public data where they didn't provide a POS file
# 
library(tidyverse)
library(optparse)
library(GenomicRanges)
# inFolder

#inFolder <- "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/WEIGHTS/PsychENCODE_eQTL"
panel <- "AMPAD_TWAS"
inFolder <- "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/WEIGHTS/AMPAD_TWAS/"
weight_folder <- "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/WEIGHTS/"
N <- 790

outFile <- paste0(inFolder, panel, ".pos")

weights <- list.files(inFolder, pattern = ".RDat", full.names = TRUE)
weights <- gsub("//", "/", weights)
weights <- gsub(weight_folder, "", weights)

# assume ensembl id for now
genes <- gsub(".wgt.RDat", "", basename(weights) )

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
df$P0 <- start(gtf_df)
df$P1 <- end(gtf_df)
df$N <- N

write_tsv(df, outFile)
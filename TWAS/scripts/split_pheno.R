# split QTL phenotype files into chunks
library(tidyverse)
library(optparse)
library(GenomicRanges)
library(rtracklayer)

option_list <- list(
    make_option(c('-i', '--input' ), help='The tensorQTL phenotype file for the dataset', default = ""),
    make_option(c('-o', '--out_prefix' ), help = "the full path to the output file for writing to", default = ""),
    make_option(c('-n', '--n_chunk' ), help = "the number of chunks to split file by", default = 100),
    make_option(c('-m', '--mode' ), help = "whether eQTLs or sQTLs are present", default = "eQTL"),
    make_option(c('-c', '--chr_type'), help = "whether chromosome should be chr1 or 1", default = "chr1"),
    make_option(c('-l', '--liftover' ), help = "whether to lift over coordinates from hg38 to hg19", default = FALSE, action = "store_true")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

input <- opt$input
out_prefix <- opt$out_prefix
mode <- opt$mode
liftover <- opt$liftover
n_chunk <- opt$n_chunk
chr_type <- opt$chr_type

liftOverCoord <- function(df, from = "hg38", to = "hg19"){
    if( from == "hg38" & to == "hg19" ){
        chain <- chain_hg38_hg19
    }else{
        stop("only hg38 -> hg19 supported currently")
    }
    coord <- data.frame(
        chr = df[[1]],
        start = df[[2]],
        end = df[[3]],
        stringsAsFactors = FALSE
        )
    #coord <- splitCoords(coord_string)
    # lift over assumes chr1 format 
    if( all(!grepl("chr", coord$chr)) ){
        coord$chr <- paste0("chr", coord$chr)
    }
    # make genomicRanges object
    # liftOver using 
    coord_gr <-  GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(coord$chr),
      ranges = IRanges::IRanges(start=as.numeric(coord$start), end=as.numeric(coord$end)),
      strand = S4Vectors::Rle(rep('*',nrow(coord)))
        )
    # lift over using rtracklayer
    # watch out for duplicate entries
    message(" * lifting over....")
    lifted_over <- rtracklayer::liftOver(coord_gr, chain )
    message(" * lifted!" )
    #return(lifted_over)

    stopifnot(length(lifted_over) == length(coord_gr)  )
    coord$chr <- seqnames(lifted_over)
    coord$start <- start(lifted_over)
    coord$end <- end(lifted_over)
    coord$chr <- gsub("chr", "", coord$chr)
    
    df[[1]] <- as.numeric(coord$chr)
    df[[2]] <- as.numeric(coord$start)
    df[[3]] <- as.numeric(coord$end)
    return(df)
}



## TODO - if the file is an sQTL, do I have to do something different here? Figure out later

df <- read_tsv(input)

# remove "chr" from chromosomes
if( chr_type != "chr1" ){
    df[[1]] <- gsub("chr", "", df[[1]])
}

if( mode == "sQTL"){
    df[[4]] <- gsub(":clu_[0-9]+_[+-]", "", df[[4]])
}

# remove tag from gene
df[[4]] <- gsub("\\.[0-9]+", "", df[[4]])

## not currently advised.
## lifting over intervals to hg19
if( liftover == TRUE){
    chain_hg38_hg19 <- import.chain("/sc/arion/projects/ad-omics/data/references/liftOver/hg38ToHg19.over.chain")
    df <- liftOverCoord(df)
}
## splitting
x <- 1:nrow(df)

message(paste0(" * splitting into ", n_chunk, " chunks") )

chunk_df <- split(df, sort(x%%n_chunk) )

# write each chunk
walk2( chunk_df, 1:n_chunk, ~{
    out_file <- paste0(out_prefix, ".pheno.", .y, ".bed")
    write_tsv(.x, out_file)
})




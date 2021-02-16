# split QTL phenotype files into chunks
library(tidyverse)
library(optparse)
option_list <- list(
    make_option(c('-i', '--input' ), help='The full path to the folder that contains the COLOC results', default = ""),
    make_option(c('-o', '--out_prefix' ), help = "the full path to the output file for writing to", default = ""),
    make_option(c('-n', '--n_chunk' ), help = "the number of chunks to split file by", default = 100),
    make_option(c('-m', '--mode' ), help = "whether eQTLs or sQTLs are present", default = "eQTL")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

input <- opt$input
out_prefix <- opt$out_prefix
mode <- opt$mode
n_chunk <- opt$n_chunk

## TODO - if the file is an sQTL, do I have to do something different here? Figure out later

df <- read_tsv(input)

x <- 1:nrow(df)

message(paste0(" * splitting into ", n_chunk, " chunks") )

chunk_df <- split(df, sort(x%%n_chunk) )

# write each chunk
walk2( chunk_df, 1:n_chunk, ~{
    out_file <- paste0(out_prefix, ".pheno.", .y, ".bed")
    write_tsv(.x, out_file)
})




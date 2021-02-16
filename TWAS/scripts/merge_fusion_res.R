# merge FUSION TWAS  output files
library(tidyverse)
library(optparse)
option_list <- list(
    make_option(c('-i', '--inFolder' ), help='The full path to the folder that contains the COLOC results', default = ""),
    make_option(c('-o', '--outFile' ), help = "the full path to the output file for writing to", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

inFolder <- opt$inFolder
outFile <- opt$outFile

message(" * writing to ", outFile)

all_files <- list.files(inFolder, pattern = "twas.txt", recursive = TRUE, full.names = TRUE )

names(all_files) <- all_files

all_res <- purrr::map_df(all_files, ~{
    read_tsv(.x)
}) %>%
    mutate(TWAS.PADJ = p.adjust(TWAS.P, method = "bonferroni") ) %>%
    arrange(TWAS.PADJ)


write_tsv(all_res, outFile)

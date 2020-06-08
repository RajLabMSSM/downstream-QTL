
library(tidyverse)
library(optparse)


option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

outFolder <- opt$outFolder

outFile <- paste0(outFolder, "/", "all_qvalue_merged.tsv")

message( " * merging qvalue results!" )
files <- list.files(path = outFolder, pattern = "qvalue.tsv",full.names = TRUE)
all <- map_df( files, read_tsv, col_types = "ccn" )

write_tsv(all, path = outFile )


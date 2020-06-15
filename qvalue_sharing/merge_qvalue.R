
library(tidyverse)
library(optparse)


option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = ""),
        make_option(c('-q', '--threshold'), help = "the significant threshold used", default = 0.05)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

outFolder <- opt$outFolder
threshold <- opt$threshold

outFile <- paste0(outFolder, "/", "all_qvalue_merged.", threshold, ".tsv")

message( " * merging qvalue results!" )
files <- list.files(path = outFolder, pattern = paste0(threshold, ".qvalue.tsv"),full.names = TRUE)
all <- map_df( files, read_tsv, col_types = "ccn" )

write_tsv(all, path = outFile )


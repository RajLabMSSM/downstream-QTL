library(optparse)
library(purrr)
library(readr)
library(dplyr)

option_list <- list(
        make_option(c('-g', '--gwas'), help='the name of the GWAS', default = ""),
        make_option(c('-q', '--qtl'), help='the name of the QTL', default = ""),
        make_option(c('-o', '--outFile'), help='the path to the output file', default = ""),
        make_option(c('-i', '--inFile'), help= "the path to the input RData file" )
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
inFile <- opt$inFile 
gwas <- opt$gwas
qtl <- opt$qtl

message(" * reading in COLOC RData from ", inFile)

load(inFile)

summary <- map_df( all_obj, ~{ map_df(.x, "df", .id = "geneid")  }, .id = "locus")

summary$GWAS <- gwas
summary$QTL <- qtl

summary <- select(summary, GWAS, QTL, everything() )

message(" * writing summary-level COLOC results to ", outFile)

write_tsv(summary, path = outFile)


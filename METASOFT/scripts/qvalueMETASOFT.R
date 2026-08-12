# qvalue meta-analysis results

# multiple testing per feature
library(optparse)
library(qvalue)
library(readr)



option_list <- list(
    make_option(c('-o', '--output'), help = "the output file", default = ""),
    make_option(c('-i', '--input'), help = "the input metasoft file", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

output <- opt$output
input <- opt$input

stopifnot(file.exists(input) )
# input - top_assoc file
# output - top_assoc with qvalue column


df <- read_tsv(input)

if("PADJ_RE2" %in% names(df) ){
    df$qvalue <- qvalue(df$PADJ_RE2)$qvalue
}


if("PADJ_RE2C" %in% names(df) ){
    df$qvalue <- qvalue(df$PADJ_RE2C)$qvalue
}

write_tsv(df, output)


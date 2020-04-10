
# Extract Top Hits
# Jack Humphrey

## Take a set of GWAS summary stats divided by chromosome
## Find all variants at a nominal P-value threshold (default 1e-5)
## Group variants into blocks (default 1 MB) and report the lowest P-value variant in that block
## Write out hits in BED format

library(dplyr)
library(readr)
library(purrr)
dir <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/processed/"
outFile <- "Nicolas_2018_hits_1e-5.tsv"
threshold <- 1e-5
blocksize <- 1e6

library(optparse)
library(readr)

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file', default = "top_hits.tsv"),
         make_option(c('-g', "--gwas_dir"), help = "The path to the directory containing the processed GWAS results by chromosome", default = "."),
        make_option(c('-t', '--threshold'), help = 'nominal significance threshold for variants', default = 1e-5),
    make_option(c('-b', '--blockSize'), help = "the genomic distance to aggregate by", default = "1e6")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
dir <- opt$gwas_dir
threshold <- as.numeric(opt$threshold)
blocksize <- as.numeric( opt$blockSize) 

# for each chr
# read in GWAS SNPs
# filter to threshold
# for the remaining SNPs
# group into 100kb blocks
# report the top hit within each block

files <- list.files(path = dir, pattern = "*_chr[0-9]*.tsv.gz$", full.names = TRUE)



top_hits <- 
    map_df(files, ~{
    print(paste0("processing file: ", .x ))
    d <- suppressMessages(read_tsv(.x))
    stopifnot( all( c("p_value", "pos") %in% names(d) ) )
    sig <- filter(d, p_value < threshold)

    sig$block <- floor(sig$pos / blocksize)

    top <- sig %>% split(.$block) %>% map_df( ~{ arrange(.x, p_value) %>% head(1) })
    return(top)
})

hits_df <- top_hits %>%
    mutate(end = pos + 1 ) %>%
    select(chr, pos, end, everything() ) %>%
    rename("#chr" = chr, start = pos )

write_tsv(hits_df, path = outFile)

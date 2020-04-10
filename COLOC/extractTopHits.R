library(dplyr)
library(readr)
library(purrr)
dir <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018"

threshold <- 1e-4
blocksize <- 1e6
# for each chr
# read in GWAS SNPs
# filter to threshold
# for the remaining SNPs
# group into 100kb blocks
# report the top hit within each block

files <- list.files(path = dir, pattern = "*_chr[0-9]*.tsv.gz$")


#d <- read_tsv("Nicolas_2018_processed_chr9.tsv.gz")


top_hits <- 
    map_df(files, ~{
    print(paste0("processing file: ", .x ))
    d <- read_tsv(.x)
    sig <- filter(d, p_value < threshold)

    sig$block <- floor(sig$pos / blocksize)

    top <- sig %>% split(.$block) %>% map_df( ~{ arrange(.x, p_value) %>% head(1) })
    return(top)
})

hits_df <- top_hits %>%
    mutate(end = pos + 1 ) %>%
    select(chr, pos, end, everything() ) %>%
    rename("#chr" = chr, start = pos )

write_tsv(outFile)

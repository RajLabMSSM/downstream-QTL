# tranpose a covariate file

library(tidyverse)
library(optparse)
option_list <- list(
    make_option(c('-i', '--input' ), help='The full path to the folder that contains the COLOC results', default = ""),
    make_option(c('-o', '--output' ), help = "the full path to the output file for writing to", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

input <- opt$input
output <- opt$output

df <- read_tsv(input)

transposed <- df %>% 
    pivot_longer(cols = -ID) %>% 
    pivot_wider(names_from = ID, values_from = value) %>%
    mutate(FID = name, IID = name) %>% 
    select(FID, IID, everything() ) %>% 
    select(-name)

# remove any spaces in covariate names
names(transposed) <- gsub(" ", "_", names(transposed) )

write_tsv(transposed, path = output)

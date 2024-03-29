library(optparse)

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file', default = ""),
        make_option(c('-i', '--inFile'), help='the path to the input file', default = ""),
        make_option(c('-m', '--mode'), help='mode - across or within?', default = "across"),
        make_option(c('-p', '--pp'), help='min PP4', default = 0.5),
        make_option(c('--debug'), help = "load all files and then save RData without running COLOC", action = "store_true", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

inFile <- opt$inFile
outFile <- opt$outFile
mode <- opt$mode
PP4 <- opt$pp
library(tidyverse)

#all_coloc <- read_tsv("../all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz")
all_coloc <- read_tsv(inFile)

all_coloc$distance <- abs(all_coloc$GWAS_pos - all_coloc$QTL_pos)
all_coloc$distance_filter <- (!is.na(all_coloc$LD) & all_coloc$LD > 0.1) | all_coloc$distance < 5e5

qtl_coloc_input <- all_coloc %>%
  filter(distance_filter == TRUE & PP.H4.abf > PP4) %>%
  select(GWAS, QTL, locus, feature) %>%
  distinct()

message(" * ", nrow(qtl_coloc_input), " phenotypes colocalize at PP4 > ", PP4)
stopifnot( nrow(qtl_coloc_input) > 0)


shared_loci <- select(qtl_coloc_input, GWAS, locus, QTL) %>%
  distinct() %>%
  group_by(GWAS, locus) %>%
  tally() %>%
  filter(n > 1)
message(" * ", nrow(shared_loci), " loci are shared between QTL types")
stopifnot( nrow(shared_loci) > 0)

qtl_coloc_input  <- filter(qtl_coloc_input, locus %in% shared_loci$locus)

# get all pairwise comparisons required
# make optional - do you want to test features within a dataset or just across?
qtl_coloc_split <- qtl_coloc_input %>%
  split(paste0(.$GWAS, ":", .$locus)) %>%
  map_df( ~{
    d <- select(.x, GWAS,locus, QTL, feature) %>% distinct() 
    qtl_pairs <- combn(unique(paste(d$QTL, d$feature) ), 2)
    locus_df <- select(d, GWAS, locus)
    pair_df <- 
      data.frame(qtl_1 = qtl_pairs[1,], qtl_2 = qtl_pairs[2,] ) %>%
      tidyr::separate(qtl_1, into = c("qtl_1", "feature_1"), sep = " ") %>%
      tidyr::separate(qtl_2, into = c("qtl_2", "feature_2"), sep = " ") %>%
      mutate(GWAS = unique(d$GWAS), locus = unique(d$locus)) %>%
      select(GWAS, locus, everything() )

  }) %>%
  distinct()

if( mode == "across"){

qtl_coloc_split <- qtl_coloc_split %>%
 mutate( within = qtl_1 == qtl_2) %>% filter(within == FALSE)
}

write_tsv(qtl_coloc_split, outFile) #"test_eQTL_edQTL_pairwise_input.tsv")



 

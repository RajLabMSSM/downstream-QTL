
# Extract Top Hits
# Jack Humphrey

## Take a set of GWAS summary stats divided by chromosome
## Find all variants at a nominal P-value threshold (default 1e-5)
## Group variants into blocks (default 1 MB) and report the lowest P-value variant in that block
## Write out hits in BED format

library(dplyr)
library(readr)
library(purrr)
library(optparse)
library(readr)

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file', default = "top_hits.tsv"),
         make_option(c('-g', "--gwas"), help = "The name of the dataset in the GWAS/QTL database", default = "VanRheenen_2021_EUR"),
        make_option(c('-t', '--threshold'), help = 'nominal significance threshold for variants', default = 1e-5),
    make_option(c('-b', '--block_size'), help = "the genomic distance to aggregate by", default = "1e6")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

pullData <- function(dataset, type = "GWAS"){
    message(Sys.time()," * selected dataset: ", dataset)
    db_path <- "/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"

    message(Sys.time()," * reading GWAS database from ", db_path)
    stopifnot( file.exists(db_path) )

    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = type, na= c("", "-","NA")))

    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )

    #stopifnot( file.exists(gwas$full_path) )

    return(gwas)
}

extractGWAS <- function(gwas, chr, refFolder = "/sc/arion/projects/ad-omics/data/references/GWAS/", force_maf = TRUE){
    # either read in config.yaml or Brian's CSV table
    gwas_path <- gwas$full_processed_path
    #gwas_path <- file.path( refFolder,  paste0(gwas$dataset, ".processed.tsv.gz" ))
    if( !file.exists(gwas_path) ){
        stop("ERROR - processed GWAS not found, make sure you ran process_GWAS.R first")
    }
    stopifnot( all(!is.na(c(gwas$full_chrom, gwas$full_chrom_type, gwas$full_pos, gwas$full_p, gwas$full_effect, gwas$full_se, gwas$full_snp, gwas$N, gwas$build) ) ))
    chrCol <- gwas$full_chrom
    posCol <- gwas$full_pos
    pvalCol <- gwas$full_p
    betaCol <- gwas$full_effect
    seCol <- gwas$full_se
    snpCol <- gwas$full_snp
    mafCol <- gwas$full_maf
    a1Col <- gwas$full_A1
    a2Col <- gwas$full_A2
    
    # make sure chromosome names in coord match GWAS chr naming
    if(gwas$full_chrom_type == "1.0" | gwas$full_chrom_type == 1){
        coord <- gsub("chr", "", chr)
    }
    if( gwas$full_chrom_type == "chr1" ){
        if(!grepl("chr", chr) ){
        coord <- paste0("chr", chr)
        }
    }

    cmd <- paste( "ml bcftools; tabix -h ", gwas_path, coord )
    message(" * running command: ", cmd)

    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    message(" * GWAS extracted!")
    stopifnot( nrow(result) > 0 )
    
    cmd <- paste( "zless ", gwas_path, " | head -1 " )
    columns <- colnames(data.table::fread( cmd = cmd ))
    col_dict <- setNames(1:length(columns), columns)

    names(result)[names(col_dict) == chrCol]  <- "chr"
    names(result)[names(col_dict) == posCol]  <- "pos"
    names(result)[names(col_dict) == pvalCol]  <- "p"
    names(result)[names(col_dict) == betaCol]  <- "effect"
    names(result)[names(col_dict) == seCol]  <- "se"
    names(result)[names(col_dict) == snpCol] <- "snp"
    names(result)[names(col_dict) == a1Col] <- "A1"
    names(result)[names(col_dict) == a2Col] <- "A2"

    return(result)
}

# 
getTopLoci <- function(gwas, block_size, threshold){
    
    sig <- gwas[gwas$p < threshold,]

    sig$block <- floor(sig$pos / block_size)

    top <- sig %>% split(.$block) %>% map_df( ~{ arrange(.x, p) %>% head(1) })
    return(top)
}

#outFile <- opt$outFile
#dir <- opt$gwas_dir
#threshold <- as.numeric(opt$threshold)
#block_size <- as.numeric( opt$block_size) 

# for each chr
# read in GWAS SNPs
# filter to threshold
# for the remaining SNPs
# group into 100kb blocks
# report the top hit within each block

dataset <- opt$gwas
block_size <- as.numeric(opt$block_size)
threshold <- as.numeric(opt$threshold)
outFile <- opt$outFile

gwas <- pullData(dataset, "GWAS")

#stop("testing")

top_hits <- map_df( 1:22, ~{
    df <- extractGWAS(gwas, chr = .x)
    top <- getTopLoci(df, block_size, threshold)
    return(top)
})


hits_df <- top_hits %>%
    dplyr::mutate(start = pos -1, end = pos ) %>%
    dplyr::mutate(locus = paste0("locus_", chr, "_", snp) ) %>%
    dplyr::select(locus, chr, start, end, snp, A1, A2, effect, se, p, pos )

# to do: match nearest gene from, errr, GENCODE?

write_tsv(hits_df, path = outFile)

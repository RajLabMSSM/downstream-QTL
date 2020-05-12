# input:

# directory containing GWAS files
# BED file of GWAS hits to test
# directory containing QTL nominal results 

# process:
# take the GWAS hit coordinates
# expand by 1MB either side
# extract the GWAS P-value distribution using tabix

# read in QTL parquet file for that chromosome
# extract SNP-Gene pairs at the coordinates
# split pairs by Gene

# for each gene:
# get together the two distributions, merge by shared position
# run COLOC to get different colocalisations
# extract top QTL SNP for locus
# record results in table - GWAS locus, Gene, top QTL SNP, COLOC probabilities
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(arrow)
library(coloc)
library(optparse)

# flank coordinates by set number of bases (default 1MB)
# work on either a coordinate string or a dataframe containing chr start and end columns
makeCoords <- function(input, flank = 1e6){
    if(all(class(input) == "character")){
        coord <- parseCoords(input)
        coord$start <- coord$start - flank
        coord$end <- coord$end + flank
        coord <- paste0(coord$chr, ":", coord$start, "-", coord$end)    
        return(coord)
    }
    if( "data.frame" %in% class(input)){
        stopifnot( all(c("chr", "start", "end") %in% names(input) ) )
        stopifnot( flank >= 0)
        coords <- paste0(input$chr, ":", input$start - flank, "-",input$end + flank) 
        return(coords)
    }
}

parseCoords <- function(coords){
    split <- as.data.frame(stringr::str_split_fixed(coords, ":|-", 3), stringsAsFactors = FALSE)
    names(split) <- c("chr", "start","end")
    split$start <- as.numeric(split$start)
    split$end <- as.numeric(split$end)
    return(split)
}

pullGWAS <- function(dataset){
    message(Sys.time()," * selected dataset: ", dataset)
    db_path <- "/sc/hydra/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"

    message(Sys.time()," * reading GWAS database from ", db_path)
    stopifnot( file.exists(db_path) )

    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = 2,na= c("", "-","NA")))

    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )

    stopifnot( all( c("full_path", "full_chrom", "full_pos") %in% names(gwas) ) )

    stopifnot( file.exists(gwas$full_path) )

    return(gwas)
}

extractLoci <- function(gwas){
    loci_path <- gwas$top_path
    # loci are stored either in comma-separated (csv), tab-separated files (txt, tsv) or in excel files (xlsx, xlsm)
    loci_file <- unlist(stringr::str_split(basename(loci_path), "\\."))
    loci_ext <- loci_file[ length(loci_file) ]
    stopifnot( loci_ext %in% c("csv", "txt", "tsv", "xlsx", "xlsm" ) )
    # read in loci file depending on file type
    if( loci_ext == "csv" ){
        loci_df <- readr::read_csv(loci_path)
    }
    if( loci_ext %in% c("txt", "tsv") ){
        loci_df <- readr::read_tsv(loci_path)
    }
    if( loci_ext %in% c("xlsx", "xlsm") ){
        stopifnot( !is.na(gwas$top_sheet) )
        stopifnot( gwas$top_sheet %in% readxl::excel_sheets(loci_path) )
        loci_df <- readxl::read_excel(loci_path, sheet = gwas$top_sheet )
    }
    #return(loci_df)
    # rename columns
    # return only core set of columns - locus, SNP, chrom, position 
    stopifnot( all(!is.na(c(gwas$top_locus, gwas$top_snp, gwas$top_chrom, gwas$top_pos) ) ) )
    loci_df$locus <- loci_df[[gwas$top_locus]] 
    loci_df$snp <- loci_df[[gwas$top_snp]]
    loci_df$chr <- loci_df[[gwas$top_chrom]]
    loci_df$pos <- loci_df[[gwas$top_pos]]
    
    loci_clean <- dplyr::select(loci_df, locus, snp, chr, pos)
    # if P-value present then include it
    if( !is.na(gwas$top_p) ){
        loci_clean$p <- loci_df[[gwas$top_p]]
    }
    # same for beta and MAF?
    return(loci_clean)
}

# load in MAF data from 1000 Genomes
# already filtered from ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/supporting/EUR.2of4intersection_allele_freq.20100804.sites.vcf.gz
# all sites with an RSID and 2+ alleles 
loadMAF <- function(path){
    message(" * loading in MAF from 1000 Genomes")
    maf <- data.table::fread(path, nThread = 4)
    names(maf) <- c("chr", "pos", "snp", "MAF")
    maf$chr <- as.character(maf$chr)
    data.table::setkey(maf, "chr")
    maf$MAF <- suppressWarnings(as.numeric(maf$MAF))
    return(maf)
}

# use 1000 Genomes MAF estimates if not present in GWAS
# match on RS ID
matchMAF <- function(data, maf){
    stopifnot( "snp" %in% names(data) )
    stopifnot( all(!is.na(data$snp) ))
    chr <- gsub("chr", "", unique(data$chr))
    stopifnot( !is.na(chr) )
    stopifnot( length(chr) == 1)
    # subset maf by chromosome
    maf <- maf[chr] 
    # match on RSID
    message(" before joining MAF: ", nrow(data), " SNPs")
    data$MAF <- maf$MAF[match(data$snp, maf$snp)]
    data <- data[ !is.na(data$MAF) ,]
    message(" after joining MAF ", nrow(data), " SNPs")
    return(data)
}

# instead of lifting over - very fiddly
# match the RSID to each SNP in a set of SNPs
matchRSID <- function(data, build){
    stopifnot( build %in% c("hg19","hg38") )
    chr <- as.character(gsub("chr", "", unique(data$chr)))
    stopifnot( !is.na(chr) )
    stopifnot( length(chr) == 1)
    # dbpsnp_hg19 and dbsnp_hg38 are data.table objects indexed by chr
    if(build == "hg19"){
        dbsnp <- as.data.frame(dbsnp_hg19[ chr ])
    }
    if(build == "hg38"){
        dbsnp <- as.data.frame(dbpsnp_hg38[chr ])
    }
    message(" before joining RSID: ", nrow(data), "SNPs")
    # here any SNP in data that has no RSID is removed
    data$snp <- dbsnp$snp[ match(data$pos, dbsnp$pos)]
    data <- data[ !is.na(data$snp), ]
    message(" after joining RSID: ", nrow(data), "SNPs")
    
    return(data)
}

# extract SNPs within coordinate range from a GWAS summary stat file
# account for different GWAS having different column naming and ordering with a GWAS_config.yaml file
# for each GWAS set the column numbers, the number of samples (N), the type of GWAS ("cc" or "quant") and the case proportion
# and whether the GWAS was hg19 or hg38
extractGWAS <- function(gwas, coord, refFolder = "/sc/hydra/projects/ad-omics/data/references/GWAS/"){
    # either read in config.yaml or Brian's CSV table
    gwas_path <- file.path( refFolder,  paste0(gwas$dataset, ".processed.tsv.gz" ))
    if( !file.exists(gwas_path) ){
        stop("ERROR - processed GWAS not found, make sure you ran process_GWAS.R first")
    }
    stopifnot( all(!is.na(c(gwas$full_chrom, gwas$full_pos, gwas$full_p, gwas$full_effect, gwas$full_se, gwas$full_snp) ) ))
    chrCol <- gwas$full_chrom
    posCol <- gwas$full_pos
    pvalCol <- gwas$full_p
    betaCol <- gwas$full_effect
    seCol <- gwas$full_se
    snpCol <- gwas$full_snp
    mafCol <- gwas$full_maf
    # assume coord is a string following chr:start-end format
    chr <- parseCoords(coord)["chr"]
    
    cmd <- paste( "ml bcftools; tabix -h ", gwas_path, coord ) 
    
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    message(" * GWAS extracted!")
    
    # get column dictionary
    cmd <- paste( "zless ", gwas_path, " | head -1 " )
    columns <- colnames(data.table::fread( cmd = cmd ))
    col_dict <- setNames(1:length(columns), columns)

    names(result)[names(col_dict) == chrCol]  <- "chr"
    names(result)[names(col_dict) == posCol]  <- "pos"
    names(result)[names(col_dict) == pvalCol]  <- "pvalues"
    names(result)[names(col_dict) == betaCol]  <- "beta"
    names(result)[names(col_dict) == seCol]  <- "varbeta"
    names(result)[names(col_dict) == snpCol] <- "snp" 
    #return(result)
    
    # deal with MAF if missing
    if( !is.na(mafCol) ){
        message(" * MAF not present - using 1000 Genomes MAF")
        result <- matchMAF(result, maf_1000g)
    }else{
    names(result)[names(col_dict) == mafCol] <- "MAF" 
    }

    # convert standard error to the variance
    result$varbeta <- result$varbeta^2
    
    result <- dplyr::select( result, snp, pvalues, beta, varbeta, MAF)
    
    result <- as.list(result)
    
    result$N <- gwas$N
    result$type <- gwas$type
    
    if( gwas$type == "cc" ){
        result$s <- gwas$prop_cases
    }
    return(result)
}

# lift over GWAS
# if GWAS mapped to hg19 then take GWAS SNPs and lift over to hg38
# column names are standardised already
liftOverGWAS <- function(gwas, chainPath = "/sc/hydra/projects/ad-omics/data/references/liftOver/hg19ToHg38.over.chain.gz"){
    # make genomicRanges object
    # liftOver using 
    gwas_gr <-  GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(gwas$chr),
      ranges = IRanges::IRanges(start=as.numeric(gwas$pos), end=as.numeric(gwas$pos)),
      strand = S4Vectors::Rle(rep('*',nrow(gwas)))
        )
    # lift over using rtracklayer
    # watch out for duplicate entries
    lifted_over <- rtracklayer::liftOver(gwas_gr, chainPath)
    stopifnot(length(lifted_over == length(gwas_gr) ) )
    gwas$chr <- seqnames(lifted_over)
    gwas$pos <- start(lifted_over)

    return(lifted_over)
}

extractQTL_parquet <- function(qtl_path, coord){
    coord_split <- parseCoords(coord)
    stopifnot(nrow(coord_split) == 1)
    chr <- unlist(coord_split["chr"])
    
    perm_file <- paste0(qtl_path, ".cis_qtl.txt.gz" )
    stopifnot(file.exists(perm_file) )
    
    # Read in permutation results
    # get out significant genes within locus coordinates
    perm_res <- readr::read_tsv(perm_file)
    perm_res <- dplyr::bind_cols(perm_res, parseCoords(perm_res$variant_id) )
    sig <- dplyr::filter(perm_res, qval < sig_level & chr ==  coord_split$chr & start >= coord_split$start & start <= coord_split$end )
    
    message( paste0( nrow(sig), " significant genes or splicing events at this locus" ) )
    if( nrow(sig) == 0 ){ return(NULL) }
    
    # read in nominal QTL associations
    parquet_file <- paste0(qtl, ".cis_qtl_pairs.", chr, ".parquet" )
    stopifnot(file.exists(parquet_file) )
    pq <- arrow::read_parquet(parquet_file)
    
    # default column numbers specify the relevant columns in the file - these may change over time 
    names(pq)[geneCol] <- "gene"
    
    pq <- dplyr::filter(pq, gene %in% sig$phenotype_id)
 

}

# extract all nominal QTL P-values overlapping the flanked GWAS hit
# split by Gene being tested
# Nominal QTL associations are stored in parquet files, one for each chromosome
# random access isn't possible (yet) so you have to read in the entire file and subset out the region of interest
extractQTL <- function(qtl, coord, sig_level = 0.05){
    # variables stored in qtl
    # if qtl type is parquet then read in parquet files
    # if qtl type is tabix then query the coordinates through tabix 
    

    # check columns
    stopifnot( all(!is.na(c(qtl$full_chrom, qtl$full_pos, qtl$full_p, qtl$full_effect, qtl$full_se, qtl$full_snp) ) ))
    chrCol <- qtl$full_chrom
    posCol <- qtl$full_pos
    pvalCol <- qtl$full_p
    betaCol <- qtl$full_effect
    seCol <- qtl$full_se
    snpCol <- qtl$full_snp
    mafCol <- qtl$full_maf
   
    # read in subset of QTLs 
    if( qtl$type == "parquet" ){
        result <- extractQTL_parquet(qtl, coord, sig_level)
    }
    if( qtl$type == "tabix" ){
        result <- extractQTL_tabix(qtl, coord, sig_level)
    }
    
    # assign column names
                 
    variant_pos <- pq[[varCol]]
    pq$chr <- chr
    pq$pos <- unlist(parseCoords( variant_pos )["start"])
    pq$snp <- paste0(chr,':', pq$pos)
    
    names(pq)[pvalCol] <- "pvalues"
    names(pq)[betaCol] <- "beta"
    names(pq)[seCol] <- "varbeta"
    pq$varbeta <- pq$varbeta^2
    names(pq)[mafCol] <- "MAF"
    # retain only associations within locus coords
    pq_subset <- 
        dplyr::select(pq, gene, snp, pos, pvalues, beta, varbeta, MAF) %>%
        dplyr::filter( pos >= coord_split$start & pos <= coord_split$end)

    # split by gene, convert to list object
    pq_split <- 
        split(pq_subset, pq_subset$gene) %>%
        purrr::map( ~{ 
            x = as.list(.x)
            x$N <- N
            x$type <- "quant"
            return(x)
        })
    return(pq_split)
}

       
# running COLOC between a GWAS locus and all eQTLs within 1MB either side
# flank coordinates
# extract GWAS nominal SNPs from region
# find GWAS hit SNP for recording
# extract QTLs within range - at set significance level
# join together and run COLOC
# output a list of objects:
## COLOC results - this should include the GWAS locus, the gene of interest, the top QTL variant and the top QTL p-value
## object - this should be a table combining the two input datasets, inner joined on "snp", to be used for plotting.
runCOLOC <- function(gwas_prefix, qtl_prefix, hit){
    coord <- parseCoords(hit)
    #coord <- makeCoords(hit, flank = 0)
    range <- makeCoords(coord, flank = 1e6)
    # extract GWAS and QTL summary for given coord
    message(" * extracting GWAS")
    g <- extractGWAS(range, gwas = gwas_prefix)
    
    # if GWAS is hg19 then lift over to hg38
    if(GWAS.build == "hg19"){
        message(" * GWAS is hg19. Lifting to hg38")
        g <- liftOverGWAS(g)
    }

    # get hit info from GWAS
    hit_snp <- paste(parseCoords(hit)[1:2],collapse = ":")
    hit_info <- as.data.frame(g, stringsAsFactors = FALSE) %>% filter(snp == hit_snp) %>% select(GWAS_SNP = snp, GWAS_P = pvalues, GWAS_beta = beta, GWAS_MAF = MAF) 
    message(" * extracting QTLs")
    q <- extractQTL(range, qtl = qtl_prefix)
    if( is.null(q) ){ return(NULL) }
    # for each gene extract top QTL SNP    
    qtl_info <- q %>% 
        map_df( ~{ 
            as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_Beta = beta, QTL_MAF = MAF)
        }) 
    message(" * running COLOC")
    coloc_res <- 
        purrr::map_df( q, ~{
            as.data.frame(t(coloc.abf(dataset1 = g, dataset2 = .x)$summary), stringsAsFactors = FALSE)
        }, .id = "gene") %>% 
        dplyr::mutate(GWAS_SNP = paste(parseCoords(hit)[1:2],collapse = ":")  ) %>% 
        dplyr::select(gene, GWAS_SNP, everything() )
    # bind coloc_res to other tables
    full_res <- full_join(qtl_info,
            full_join(hit_info, coloc_res, by = "GWAS_SNP"), 
                by = "gene") %>% 
        select( gene, starts_with("GWAS"), starts_with("QTL"), everything() )
    return(full_res) 
}

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file', default = ""),
        make_option(c('--hits'), help= "the path to a file containing GWAS hits" ),
        make_option(c('-p', '--gwas_prefix'), help = "the prefix of the GWAS files"),
        make_option(c('-q','--qtl_prefix'), help = "the directory containing QTL nominal results with the prefix of the dataset; eg /results/Brain_expression/peer30/Brain_expression_peer30", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
hits_file <- opt$hits
qtl_prefix <- opt$qtl_prefix
gwas_prefix <- opt$gwas_prefix


#gwas_prefix <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/processed/Nicolas_2018_processed_"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_expression/peer30/LumbarSpinalCord_expression_peer30"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_splicing/peer20/LumbarSpinalCord_splicing_peer20"

#hits_file <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/COLOC/Nicolas_2018/Nicolas_2018_hits_1e-7.tsv"
#outFile <- "test_coloc_results.tsv"

#options(echo = TRUE)

main <- function(){


maf_1000g <- loadMAF("/sc/hydra/projects/ad-omics/data/references/1000G/1000G_EUR_MAF.bed.gz")

dbsnp_hg19_path <- "/sc/hydra/projects/ad-omics/data/references/hg19_reference/dbSNP/hg19_common_snps_dbSNP_153.bed.gz"
dbsnp_hg38_path <- "/sc/hydra/projects/ad-omics/data/references/hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz"

dbsnp_hg19 <- loadDBSNP(dbsnp_hg19_path)
dbsnp_hg38 <- loadDBSNP(dbsnp_hg38_path)




hits <- read_tsv(hits_file)
hit_coords <- makeCoords(hits, flank = 0)
# for testing
#hit_coords <- hit_coords[1]

all_res <- purrr::map_df(1:length(hit_coords), ~{
    message(hit_coords[.x])
    res <- runCOLOC(gwas_prefix, qtl_prefix, hit = hit_coords[.x])
    if(is.null(res) ){return(NULL) }
    return(res)
})


readr::write_tsv(all_res, path = outFile)
}

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
suppressPackageStartupMessages(library(tidyverse))

suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(optparse))

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

options(future.globals.maxSize = 8000 * 1024^2)

#library(arrow)

message_ <- function(text){
    message(Sys.time(), " * ", text)   
}
## verbose messaging
message_verbose <- function(input){
    if(verbose){
        if( "data.frame" %in% class(input)){
            print(input)
        }
        if( class(input) == "character"){
            message(Sys.time(), "     * ", input)
        }
    }
}
# flank coordinates by set number of bases (default 1MB)
# work on either a coordinate string or a dataframe containing chr start and end columns
joinCoords <- function(input, flank = 0){
    if(all(class(input) == "character")){
        coord <- splitCoords(input)
        coord$start <- coord$start - flank
        coord$end <- coord$end + flank
        if( coord$start < 0){coord$start <- 0}
        coord <- paste0(coord$chr, ":", coord$start, "-", coord$end)    
        return(coord)
    }
    if( "data.frame" %in% class(input)){
        stopifnot( all(c("chr", "start", "end") %in% names(input) ) | all(c("chr","pos") %in% names(input) )  )
        stopifnot( flank >= 0)
        if( all( c("chr","start","end") %in% names(input) ) ){
            # deal with flanks that go negative
            starts <- input$start - flank
            starts[starts < 0 ] <- 1
            coords <- paste0(input$chr, ":", starts, "-",input$end + flank) 
        }
        if( all( c("chr", "pos") %in% names(input) ) ){
            starts <- input$pos - (flank + 1)
            starts[starts < 0 ] <- 1
            coords <- paste0(input$chr, ":", starts, "-", input$pos + flank)
        }
        return(coords)
    }
}

splitCoords <- function(coords){
    split <- as.data.frame(stringr::str_split_fixed(coords, ":|-", 3), stringsAsFactors = FALSE)
    names(split) <- c("chr", "start","end")
    split$start <- as.numeric(split$start)
    split$end <- as.numeric(split$end)
    return(split)
}

pullData <- function(dataset, type = "GWAS"){
    message_verbose(paste0("selected dataset: ", dataset))
    db_path <- "/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"

    message_verbose(paste0("reading GWAS database from ", db_path))
    stopifnot( file.exists(db_path) )

    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = type, na= c("", "-","NA")))
    
    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )


    #stopifnot( file.exists(gwas$full_path) )

    return(gwas)
}

extractLoci <- function(gwas){
    loci_path <- gwas$top_filtered_path
    
    # stupid exception for Daner_2020 - filtered loci are totally different file right now
    if(gwas$dataset == "Daner_2020"){
        loci_path <- gwas$top_path
    }
    # loci are stored either in comma-separated (csv), tab-separated files (txt, tsv) or in excel files (xlsx, xlsm)
    loci_file <- unlist(stringr::str_split(basename(loci_path), "\\."))
    loci_ext <- loci_file[ length(loci_file) ]
    stopifnot( loci_ext %in% c("csv", "txt", "tsv", "xlsx", "xlsm", "xls" ) )
    # read in loci file depending on file type
    if( loci_ext == "csv" ){
        loci_df <- readr::read_csv(loci_path)
    }
    if( loci_ext %in% c("txt", "tsv") ){
        loci_df <- readr::read_tsv(loci_path)
    }
    if( loci_ext %in% c("xlsx", "xlsm", "xls") ){
        stopifnot( !is.na(gwas$top_sheet) )
        stopifnot( gwas$top_sheet %in% readxl::excel_sheets(loci_path) )
        loci_df <- readxl::read_excel(loci_path, sheet = gwas$top_sheet )
    }
    #return(loci_df)
    # rename columns
    # return only core set of columns - locus, SNP, chrom, position 
    stopifnot( all(!is.na(c(gwas$top_locus, gwas$top_snp, gwas$top_chrom, gwas$top_pos) ) ) )
    
    stopifnot( all( c(gwas$top_locus, gwas$top_snp, gwas$top_chrom) %in% colnames(loci_df) ) )
    loci_df$locus <- loci_df[[gwas$top_locus]] 
    loci_df$snp <- loci_df[[gwas$top_snp]]
    loci_df$chr <- loci_df[[gwas$top_chrom]]
    loci_df$pos <- as.numeric(loci_df[[gwas$top_pos]])

    # remove "chr" from chr name if present
    loci_df$chr <- gsub("chr", "", loci_df$chr) 
    # Ripke has 23 (chrX) in loci - remove 
    loci_df <- dplyr::filter(loci_df, chr %in% c(1:22) | chr %in% paste0("chr", 1:22) ) 
 
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
    #message_("loading in MAF from 1000 Genomes")
    maf <- data.table::fread(path, nThread = 4)
    names(maf)[1:4] <- c("chr", "pos", "snp", "MAF")
    maf$chr <- as.character(maf$chr)
    data.table::setkey(maf, "chr")
    maf$MAF <- suppressWarnings(as.numeric(maf$MAF))
    return(maf)
}

# use 1000 Genomes MAF estimates if not present in GWAS
# match on RS ID
matchMAF <- function(data, maf){
    stopifnot( "snp" %in% names(data) )
    stopifnot( "chr" %in% names(data) )
    #stopifnot( all(!is.na(data$snp) ))
    chr <- gsub("chr", "", unique(data$chr))
    stopifnot( !is.na(chr) )
    stopifnot( length(chr) == 1)
    # subset maf by chromosome
    maf <- maf[chr] 
    # match on RSID
    message_verbose(paste0("before joining MAF: ", nrow(data), " SNPs"))
    data$MAF <- maf$MAF[match(data$snp, maf$snp)]
    matched <- data[ !is.na(data$MAF) & data$MAF >0 & data$MAF < 1,]
    message_verbose(paste0("after joining MAF: ", nrow(matched), " SNPs"))
    if(nrow(matched) == 0 ){
        #message_verbose(head(data) )
        message_verbose("no SNPs match on RSID. Inspect data")
        return(NULL)
    }
    return(matched)
}

# 
# extract SNPs within coordinate range from a GWAS summary stat file
# account for different GWAS having different column naming and ordering with a GWAS_config.yaml file
# for each GWAS set the column numbers, the number of samples (N), the type of GWAS ("cc" or "quant") and the case proportion
# and whether the GWAS was hg19 or hg38
extractGWAS <- function(gwas, coord, refFolder = "/sc/arion/projects/ad-omics/data/references/GWAS/", force_maf = TRUE){
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
    # assume coord is a string following chr:start-end format
    chr <- splitCoords(coord)["chr"]
    # make sure chromosome names in coord match GWAS chr naming
    if(gwas$full_chrom_type == "1.0"){
        coord <- gsub("chr", "", coord)
    }
    if( gwas$full_chrom_type == "chr1" ){
        if(!grepl("chr", coord) ){
        coord <- paste0("chr", coord)
        }
    }
     
    cmd <- paste( "ml bcftools/1.9; tabix -h ", gwas_path, coord ) 
    message_verbose(paste0("running command: ", cmd) )
 
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    message_verbose("GWAS extracted!")
    stopifnot( nrow(result) > 0 )
    
     
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
    names(result)[names(col_dict) == a1Col] <- "A1"
    names(result)[names(col_dict) == a2Col] <- "A2"
    #return(result)
    
    # deal with MAF if missing
        if( is.na(mafCol) | force_maf == TRUE ){
            message_verbose("MAF not present - using 1000 Genomes MAF")
            # how to get this in to the function? it can't find maf_1000g
            result <- matchMAF(result, maf = maf_1000g)
        }else{
            message_verbose("using supplied MAF")
            names(result)[names(col_dict) == mafCol] <- "MAF" 
        }
    # convert standard error to the variance
    # se is standard error - not standard deviation!
    result$varbeta <- result$varbeta^2
    #print(head(result))
    #save.image("debug.RData") 
    result <- dplyr::select( result, snp, pvalues, beta, varbeta, MAF, chr, pos, A1, A2)
    
    result <- as.list(result)
    
    result$N <- as.numeric(gwas$N)
    result$type <- gwas$type
    
    if( gwas$type == "cc" ){
        if( !is.na(gwas$prop_cases)){
        result$s <- as.numeric(gwas$prop_cases)
        }else{
            result$s <- 0.5
        }
    }
    return(result)
}

# lift over coordinate string
# if GWAS mapped to hg19 then take GWAS locus coordinate and lift over to hg38
liftOverCoord <- function(coord_string, from = "hg19", to = "hg38"){
    if( from == "hg19" & to == "hg38" ){
        chain <- chain_hg19_hg38
    }else{
        stop("only hg19 -> hg38 supported currently")
    }
    coord <- splitCoords(coord_string)
    # lift over assumes chr1 format 
    if( !grepl("chr", coord$chr) ){
        coord$chr <- paste0("chr", coord$chr)
    } 
    # make genomicRanges object
    # liftOver using 
    coord_gr <-  GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(coord$chr),
      ranges = IRanges::IRanges(start=as.numeric(coord$start), end=as.numeric(coord$end)),
      strand = S4Vectors::Rle(rep('*',nrow(coord)))
        )
    # lift over using rtracklayer
    # watch out for duplicate entries
    message_verbose("lifting over....")
    lifted_over <- rtracklayer::liftOver(coord_gr, chain )
    message_verbose("lifted!" )
    #return(lifted_over)
    
    stopifnot(length(lifted_over) == length(coord_gr)  )
    coord$chr <- seqnames(lifted_over)
    coord$start <- start(lifted_over)
    coord$end <- end(lifted_over)
    coord$chr <- gsub("chr", "", coord$chr)
    
    lifted_string <- joinCoords(coord, flank = 0)
    return(lifted_string)
}

# if QTL nominal stats stored in parquet format, split by chromosome
extractQTL_parquet <- function(qtl, coord, sig_level = 0.05){
    coord_split <- splitCoords(coord)
    stopifnot(nrow(coord_split) == 1)
    chr <- unlist(coord_split["chr"])
    
    if( !grepl("chr", chr) ){ chr <- paste0("chr", chr) }    
 
    perm_file <- qtl$top_path
    stopifnot(file.exists(perm_file) )
    
    # Read in permutation results
    # get out significant genes within locus coordinates
    perm_res <- readr::read_tsv(perm_file)
    perm_res <- dplyr::bind_cols(perm_res, splitCoords(perm_res$variant_id) )
    sig <- dplyr::filter(perm_res, qval < sig_level & chr ==  coord_split$chr & start >= coord_split$start & start <= coord_split$end )
    
    message_verbose( paste0( nrow(sig), " significant genes or splicing events at this locus" ) )
    if( nrow(sig) == 0 ){ return(NULL) }
    
    # read in nominal QTL associations
    parquet_file <- paste0(qtl$full_path, ".cis_qtl_pairs.", chr, ".parquet" )
    stopifnot(file.exists(parquet_file) )
    pq <- arrow::read_parquet(parquet_file)
    
    # default column numbers specify the relevant columns in the file - these may change over time 
    names(pq)[geneCol] <- "gene"
    
    pq <- dplyr::filter(pq, gene %in% sig$phenotype_id)
 
    return(pq)
}

# if QTL nominal stats stored in single tabixed file
extractQTL_tabix <- function(qtl, coord){
    if( qtl$full_chrom_type == "chr1" ){
        if( !grepl("chr", coord) ){ 
            coord <- paste0("chr", coord)
        }
    }
    stopifnot( file.exists(qtl$full_path) )
    cmd <- paste( "ml bcftools/1.9; tabix -h ", qtl$full_path, coord )
    message_verbose(paste0("running command: ", cmd))
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    # deal with regions of no QTL association
    if( nrow(result) == 0){ 
        return(NULL) 
    }
    #stopifnot( nrow(result) > 0 )
    message_verbose("QTL extracted!")
    return(result)
}


# extract all nominal QTL P-values overlapping the flanked GWAS hit
# split by Gene being tested
# sig_level doesn't do anything yet
extractQTL <- function(qtl, coord, sig_level = 0.05, force_maf = FALSE, targets = NULL){
    # variables stored in qtl
    # if qtl type is parquet then read in parquet files
    # if qtl type is tabix then query the coordinates through tabix 
    

    # check columns
    stopifnot( all(!is.na(c(qtl$full_snp, qtl$full_pheno, qtl$full_p, qtl$N, qtl$build) ) ))
    pvalCol <- qtl$full_p
    betaCol <- qtl$full_effect # if not present will be set to NA
    phenoCol <- qtl$full_pheno
    seCol <- qtl$full_se
    snpCol <- qtl$full_snp
    mafCol <- qtl$full_maf
    chrCol <- qtl$full_chrom
    posCol <- qtl$full_pos
    # read in subset of QTLs 
    if( qtl$full_file_type == "parquet" ){
        result <- extractQTL_parquet(qtl, coord, sig_level)
    }
    if( qtl$full_file_type == "tabix" ){
        result <- extractQTL_tabix(qtl, coord)
    }
    # if no SNPs found in region return NULL value
    if(is.null(result) ){ return(NULL) }
        
    # assign column names
    cmd <- paste( "zless ", qtl$full_path, " | head -1 " )
    columns <- colnames(data.table::fread( cmd = cmd ))
    col_dict <- setNames(1:length(columns), columns)
    
    #print(col_dict)
    #stopifnot( ncol(result) == length(col_dict) )
    col_dict <- col_dict[1:ncol(result)] # hacky workaround for now
    names(result)[names(col_dict) == phenoCol] <- "gene"
    names(result)[names(col_dict) == mafCol]   <- "MAF"
    names(result)[names(col_dict) == pvalCol]  <- "pvalues"
    if( !is.na(betaCol) ){
        names(result)[names(col_dict) == betaCol]  <- "beta"
    }
    if( !is.na(seCol) ){
        names(result)[names(col_dict) == seCol]    <- "varbeta"
        # don't forget to square the standard error to get the variance             
        result$varbeta <- result$varbeta^2
        message_verbose(" remove associations with 0 variance")
        message_verbose(paste0(" before: ", nrow(result), " associations"))
        # remove rows with variance of 0 - this will corrupt COLOC
        result <- dplyr::filter(result, varbeta != 0) 
        message_verbose(paste0(" after: ", nrow(result), " associations"))
    }
    names(result)[names(col_dict) == snpCol]   <- "snp"
    # Young microglia use log10P - convert
    if( pvalCol == "log10_p"){
        result$pvalues <- 10^result$pvalues
    }

    # deal with MAF - meta-analysis outputs won't have it
    # this requires "chr" to be present in QTL data
    #if( debug == TRUE){ result$MAF <- NA }else{ # why do we do this? breaks it when debug is TRUE
    if( is.na(mafCol) | force_maf == TRUE ){
        message_verbose("MAF not present - using 1000 Genomes MAF")
        #stopifnot( "chr" %in% names(col_dict) )
        names(result)[names(col_dict) == "chr"] <- "chr"
        result <- matchMAF(result, maf = maf_1000g)
        if( all(is.na(result$MAF) ) | is.null(result) ){ return(NULL) }
        #print(head(result) )
    }else{
        message_verbose("using supplied MAF")
        names(result)[names(col_dict) == mafCol] <- "MAF"
        #print(col_dict)
        #print(head(result) )
    }
    #}
    # deal with edge cases - 1 gene and the gene id is NA
    if( all( is.na(result$gene) ) ){
        message_verbose("no genes in QTL result")
        return(NULL)
    }else{
        result <- dplyr::filter(result, !is.na(gene) )
    }
    # if MAF needs matching then chrCol will be "chr" - change it back to "QTL_chr"
    #print( names(col_dict) )
    #print( names(col_dict) == chrCol )
    names(result)[which(names(col_dict) == chrCol) ]   <- "QTL_chr"
    names(result)[which(names(col_dict) == posCol) ]   <- "QTL_pos"

    #print(head(result) ) 
    
    # retain only associations within locus coords
    if( !is.na(betaCol) & !is.na(seCol) ){
        res_subset <- dplyr::select(result, gene, snp, pvalues, beta, varbeta, MAF, QTL_chr, QTL_pos)
    }else{
        res_subset <- dplyr::select(result, gene, snp, pvalues, MAF, QTL_chr, QTL_pos)
    }
    # jan 2025: another edge case - NA pvalues
    message_verbose("removing NA P-values")
    res_subset <- filter(res_subset, !is.na(pvalues) )
    message_verbose(paste0("after: ", nrow(res_subset) ))
    # edge case - multiallelic SNPs lead to two entries per SNP-gene - remove
    message_verbose("removing multi-allelic SNPs")
    res_subset <- dplyr::distinct(res_subset)
    message_verbose(paste0("after: ", nrow(res_subset) ))
    #print(head(res_subset) ) 
    message_verbose("passed this point")
    #save(list = ls(all.names = TRUE), file = "image.RData", envir = 
    #environment())
    #dplyr::filter( pos >= coord_split$start & pos <= coord_split$end)
      
    if( !is.null(targets) ){
        res_subset <- dplyr::filter(res_subset, gene == targets)
        message_verbose(paste0(length(unique(res_subset$gene)), " target features remain" ) ) 
        message_verbose(paste0(nrow(res_subset), " associations remain"))
        if( nrow(res_subset) == 0 ){
            return(NULL)
        }
    }
    # split by gene, convert to list object
    res_split <- 
        split(res_subset, res_subset$gene) %>%
        purrr::map( ~{ 
            x = as.list(.x)
            x$N <- as.numeric(qtl$N)
            x$type <- "quant"
            return(x)
        })
    return(res_split)
}

# for each gene extract top QTL SNP   
getQTLinfo <- function(q, hit_info){  
    if( all( c("beta", "varbeta") %in% names(q[[1]]) ) ){
        qtl_info <- q %>% 
        map_df( ~{ 
            as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_Beta = beta, QTL_SE = varbeta, QTL_MAF = MAF, QTL_chr, QTL_pos) %>%
            mutate( QTL_SE = sqrt(QTL_SE) ) # get SE back from Variance
        })
        # get out beta, se and P for GWAS SNP 
        gwas_info <- q %>%
         map_df( ~{
            df <- as.data.frame(.x, stringsAsFactors=FALSE) %>%
                filter(snp == hit_info$GWAS_SNP )
            if( nrow(df) == 1){
                df %>%
                select(GWAS_SNP_Beta = beta,
                       GWAS_SNP_SE = varbeta,
                       GWAS_SNP_P = pvalues ) %>%
                mutate(GWAS_SNP_SE = sqrt(GWAS_SNP_SE) )
            }else{
            tibble(GWAS_SNP_Beta = NA, GWAS_SNP_SE = NA, GWAS_SNP_P = NA)
            }
        }) 
        qtl_info <- bind_cols(qtl_info, gwas_info)
        message_verbose(head(qtl_info) )
    }else{
        qtl_info <- q %>% 
        map_df( ~{ 
            as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_MAF = MAF, QTL_chr, QTL_pos)
        }) 

    }
    return(qtl_info)
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
## sig.level now filters out QTLs with P > sig.level
## bare - if TRUE then only return COLOC summary, not full object (saves memory)
runCOLOC <- function(gwas, qtl, qtl2 = NULL, hit, sig.level = NULL, target.file = NULL, bare = FALSE){
    message_(paste0("Analysing locus: ", hit$locus ))
    # hit is a dataframe containing "snp", "chr", "pos", "locus"
    
    hit_coord <- joinCoords(hit, flank = 0)
    hit_range <- joinCoords(hit, flank = 1e6)
 
    # extract GWAS and QTL summary for given coord
    if( is.null(qtl2) ){
        message_verbose("extracting GWAS")
        g <- extractGWAS(gwas, hit_range)
    }
    # get hit info from GWAS
    hit_snp <- hit$snp
  
    # problem when no p column in hit!
    if(is.null(gwas) ){
        hit_info <- hit %>%
            select(locus) %>% mutate(GWAS_SNP = NA, GWAS_P = NA)
        qtl_range <- hit_range
    }else{
        hit_info <- hit %>%
            select(locus, GWAS_SNP = snp, GWAS_P = p)
        # if GWAS is hg19 then lift over to hg38
        if(gwas$build != qtl$build){
            qtl_coord <- liftOverCoord(hit_coord, from = gwas$build, to = qtl$build)
        }else{
            qtl_coord <- hit_coord
        }
        # record the position of the GWAS top SNP
        hit_info$GWAS_chr <- splitCoords(qtl_coord)$chr
        hit_info$GWAS_pos <- splitCoords(qtl_coord)$end
        qtl_range <- joinCoords(qtl_coord, flank = 1e6)
    }
    message_verbose("extracting QTLs")
    q <- extractQTL(qtl, qtl_range, targets = target.file[1])
    #print(names(q[[1]]))
    if( is.null(q) ){ return(NULL) }
    qtl_info <- getQTLinfo(q, hit_info)
     # if sig filtering requested
    if( !is.null(sig.level ) ){
        qtl_info <- filter(qtl_info, QTL_P < sig.level)
        q <- q[ qtl_info$gene ] 
         message_verbose(paste0(length(q), " features in QTL"))
    }

    if( !is.null(qtl2) ){
        # if second QTL dataset requested
        q2 <- extractQTL(qtl2, qtl_range, targets = target.file[2])
        if( is.null(q2) ){ return(NULL) }
        # get QTL betas for lead SNP in QTL 1
        hit_info_qtl1 <- data.frame(locus = hit_info$locus, GWAS_SNP = qtl_info$QTL_SNP, GWAS_P = qtl_info$QTL_P)
        qtl2_info <- getQTLinfo(q2, hit_info_qtl1)
        # get QTL betas for lead SNP in QTL 2
        hit_info_qtl2 <- data.frame(locus = hit_info$locus, GWAS_SNP = qtl2_info$QTL_SNP, GWAS_P = qtl2_info$QTL_P)
        qtl_info <- getQTLinfo(q, hit_info_qtl2)
        if( !is.null(sig.level ) ){
            qtl2_info <- filter(qtl2_info, QTL_P < sig.level)
            q2 <- q2[ qtl2_info$gene ]
        }
        message_verbose( paste0(length(q2), " features in QTL 2"))
    }
    # actually run COLOC
    message_("running COLOC")
    # return the results table and an object containing g and q
    # GWAS-QTL COLOC
    if( is.null(qtl2) ){
    message_("GWAS-QTL COLOC")
    coloc_res <- 
        purrr::map( q, ~{
            if( length(intersect(g$snp, .x$snp) ) == 0 ){ 
                message("Hold up! GWAS and QTL have no SNPs in common")
                return(NULL) 
            }
            coloc_object <- coloc.abf(dataset1 = g, dataset2 = .x)
            # add in g and q to coloc object
            gq <- dplyr::inner_join(as.data.frame(g),as.data.frame(.x), by = "snp", suffix = c(".gwas", ".qtl") )
            coloc_object$results <- dplyr::inner_join(coloc_object$results, gq, by = "snp" )
            coloc_df <- as.data.frame(t(coloc_object$summary), stringsAsFactors = FALSE)
            return( list(df = coloc_df, object = coloc_object) )
        })
    }else{
    message_("QTL-QTL COLOC")
    # remove features with lead SNPs > 1e-5
    # test all pairs of QTL features in the locus
    qtl_combos <- expand.grid(1:length(q), 1:length(q2) )
    message_(paste0("testing ", nrow(qtl_combos), " pairs of QTLs"))
    # qtl_info is matched on to COLOC results using gene names - but for QTL-QTL there are two genes, so instead use the index of the pair
    qtl_info <- 
        full_join( 
            qtl_info[qtl_combos$Var1,] %>% mutate(index = 1:nrow(qtl_combos) ), 
            qtl2_info[qtl_combos$Var2,] %>% mutate(index = 1:nrow(qtl_combos) ), 
            by = "index", 
            suffix = c(".qtl1", ".qtl2") ) %>%
            dplyr::mutate(gene = as.character(index)) %>%
            select(-index)

    coloc_res <- 
        purrr::map( 1:nrow(qtl_combos), ~{
            q_1 = q[[ qtl_combos[.x, "Var1"] ]]
            q_2 = q2[[ qtl_combos[.x, "Var2"] ]]
            message_( paste0( .x, ": ", unique(q_1$gene), " vs ",  unique(q_2$gene) ) )
            if( length(intersect(q_1$snp, q_2$snp) ) == 0 ){
                message_("Hold up! GWAS and QTL have no SNPs in common")
                return(NULL)
            }
            coloc_object <- coloc.abf(dataset1 = q_1, dataset2 = q_2)
            qq <- dplyr::inner_join(as.data.frame(q_1),as.data.frame(q_2), by = "snp", suffix = c(".qtl1", ".qtl2") )
            coloc_object$results <- dplyr::inner_join(coloc_object$results, qq, by = "snp" )
            coloc_df <- as.data.frame(t(coloc_object$summary), stringsAsFactors = FALSE)
            return( list(df = coloc_df, object = coloc_object) )
        })
    }
    if( is.null(unlist(coloc_res)) ){ return(NULL) }

    message_("COLOC finished")
    coloc_df <- map_df(coloc_res, "df", .id = "gene") %>% 
        dplyr::mutate(locus = hit$locus ) %>%
        dplyr::select(locus, gene, everything() )
    message_("COLOC res created")
    # bind coloc_res to other tables
    full_res <- full_join(qtl_info, coloc_df, by = "gene")
    full_res <- full_join(hit_info, full_res, by = "locus")
    # remove rows with no coloc run
    full_res <- filter(full_res, !is.na(PP.H4.abf) )
    if(bare == TRUE){
        return(list(full_res = full_res) )
    }else{
        return(list(full_res = full_res, coloc_object = coloc_res) ) 
    }
}

# for each COLOC result per locus
# calculate pairwise LD between GWAS SNP and QTL top SNPs
calc_LD <- function( coloc_res ){
    require(LDlinkR)
    # get LDlink token from .Renviron (only Jack has this, for access go to https://ldlink.nci.nih.gov/?tab=apiaccess)
    token <- Sys.getenv("LDLINK_TOKEN")
    if( token == ""){ 
        warning(Sys.time(), " * LDlink token not found!" )
        return(NA)
    }
    snps <- unique( c( unique(coloc_res$GWAS_SNP), unique( coloc_res$QTL_SNP) ) )
    
    # remove NA values
    snps <- snps[!is.na(snps)]
    
    # get pairwise LD matrix
    ld_matrix <-  LDmatrix( snps, pop = "CEU", r2d = "r2", token = Sys.getenv("LDLINK_TOKEN") )
    
    #stopifnot( nrow(ld_matrix) < length(snps) )

    ld_matrix <- tibble::column_to_rownames(ld_matrix, var = "RS_number")
    
    coloc_res$LD_CEU_R2 <- ld_matrix[ coloc_res$QTL_SNP, coloc_res$GWAS_SNP][,1]
    # wait 5 seconds before returning - makes sure API queries are spread out
    #Sys.sleep(time = 5)
    return(coloc_res) 
}


option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = "."),
        make_option(c('--gwas', '-g'), help= "the dataset ID for a GWAS in the GWAS/QTL database", default = NULL ),
        make_option(c('--qtl', '-q'), help = "the dataset ID for a QTL dataset in the GWAS/QTL database"),
        make_option(c('--qtl2', '-r'), help = "a second QTL dataset - using this will trigger QTL-QTL COLOC", default = NULL),
        make_option(c('--sig', '-s'), help = "the minimum p-value a QTL should have for testing", default = NULL),
        make_option(c('--loci', '-l'), help = "text file of loci names", default = NULL),
        make_option(c('--targets', '-t'), help = "TSV file of features, loci and GWAS for targeted analysis", default = NULL),
        make_option(c('--debug'), help = "load all files and then save RData without running COLOC", action = "store_true", default = FALSE),
        make_option(c('--lowmem', '-b'), help = "do not save COLOC objects, just the summaries", action = "store_true", default = FALSE),
        make_option(c('--threads', '-c'), help = "how many threads to use", default = 1),
        make_option(c('--MAF', '-m'), help = "path to custom MAF file", default = NULL),
        make_option(c('--verbose', '-v'), help = "if selected then output more information", action = "store_true", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

outFolder <- opt$outFolder
gwas_dataset <- opt$gwas
qtl_dataset <- opt$qtl
qtl_dataset_2 <- opt$qtl2
debug <- opt$debug
sig <- opt$sig
loci <- opt$loci
target_file <- opt$targets
low_mem <- opt$lowmem
threads <- opt$threads
verbose <- opt$verbose
MAF <- opt$MAF

# default MAF is 1000g phase 3v5
# but users can provide their own
# must have columns chr, pos, snp, MAF
if( is.null(MAF) ){
    #maf_1000gp1 <- "/sc/arion/projects/ad-omics/data/references/1KGP1/1000G_EUR_MAF.bed.gz"
    MAF <- "/sc/arion/projects/ad-omics/data/references/1KGPp3v5/EUR_MAF/EUR.all.phase3_MAF.bed.gz"
}

if(verbose){print(opt)}

# load in MAF table
if( !exists("maf_1000g") ){
  message(paste0("loading MAF from ", MAF ) )
  maf_1000g <- loadMAF(MAF)
  # load in liftover chain
  chain_hg19_hg38 <- import.chain("/sc/arion/projects/ad-omics/data/references/liftOver/hg19ToHg38.over.chain")    
}

main <- function(){
    message_("starting COLOC pipeline")
    if( debug == TRUE){ save.image("debug.RData") }
    # regular GWAS-QTL COLOC
    if( is.null(target_file) ){ 
        message_("GWAS - QTL COLOC")
        gwas <- pullData(gwas_dataset, "GWAS")
        qtl <- pullData(qtl_dataset, "QTL")
    
        if( !is.null(qtl_dataset_2) ){
            qtl2 <- pullData(qtl_dataset_2, "QTL")
        }else{ qtl2 <- NULL}
        top_loci <- extractLoci(gwas)
        message_(paste0(" testing ", nrow(top_loci), " loci ") )

        all_coloc <- purrr::map(
        1:nrow(top_loci), ~{
            res <- runCOLOC(gwas, qtl, qtl2, hit = top_loci[.x,], sig.level = sig)
            if(is.null(res) ){return(NULL) }
            return(res)
        })
    }
    # QTL-QTL COLOC with pairwise targets
    if(!is.null(target_file) ){
        message_("QTL-QTL COLOC")
        targets <- read_tsv(target_file)
        # for TESTING
        #targets <- head(targets, 3)
        message_( paste0("testing ", nrow(targets), " pairs of QTLs" ) )
        # doesn't work properly - furrr loads the whole MAF object into memory for each parallel run, 100MB per run which is crazy.
        # for parallel to work I would have to transition to querying the MAF database each time like we do with the sumstats. 
        if(threads > 1){
            suppressPackageStartupMessages(library(furrr))
            plan(multisession, workers = threads)
            map <- future_map
        }
        all_coloc <- map(
        1:nrow(targets), ~{
            message_verbose( paste0("Pair ", .x , " of ", nrow(targets) ) )
            pair_df <- targets[.x,]
            if( "GWAS" %in% names(pair_df) ){
                gwas <- pullData(pair_df$GWAS, "GWAS") 
                locus <- extractLoci(gwas) %>% filter(locus == pair_df$locus)
            }else{
                message_verbose("no GWAS specified - using chr and pos from targets file")
                locus <- pair_df
                gwas <- NULL
            }
            qtl1 <- pullData(pair_df$qtl_1, "QTL")
            qtl2 <- pullData(pair_df$qtl_2, "QTL")
            features <- c(pair_df$feature_1, pair_df$feature_2)
            res <- runCOLOC(gwas, qtl1, qtl2, hit = locus, sig.level = sig, target.file = features, bare = low_mem)
            gc()
            if(is.null(res) ){return(NULL) }
            return(res)
        }
        )
    }
    
    all_res <- purrr::map_df(all_coloc, "full_res")
    if(!is.null(target_file) ){
        outFile <- file.path(outFolder, gsub(".tsv", "_COLOC.tsv", basename(target_file) ) )
        # here if any of the COLOCs did not proceed then the target df will be larger than the all_res results        

        #all_res <- dplyr::bind_cols(targets, all_res %>% select(-locus) )
    }else{
        outFile <- paste0(outFolder,"/", qtl_dataset, "_", gwas_dataset, "_COLOC.tsv")
        names(all_obj) <- top_loci$locus
    }
    dir.create(outFolder,  showWarnings = FALSE)
    all_res <- dplyr::arrange(all_res, desc(PP.H4.abf)  )
    message("finished COLOC pipeline")
    message_verbose(paste0("writing output to ", outFile))
    readr::write_tsv(all_res, file = outFile)
    # save COLOC objects for plotting
    if( !low_mem ){
        all_obj <- purrr::map(all_coloc, "coloc_object")
        if( is.null(targets) ){
            names(all_obj) <- top_loci$locus
        }
        outData <- gsub("tsv.gz", "RData", outFile)
        outData <- gsub("tsv", "RData", outFile)
        save(all_obj, file = outData)
    }
}

main()



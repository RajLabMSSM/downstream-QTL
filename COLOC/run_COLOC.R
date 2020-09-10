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
library(coloc)
library(optparse)
#library(arrow)

# flank coordinates by set number of bases (default 1MB)
# work on either a coordinate string or a dataframe containing chr start and end columns
joinCoords <- function(input, flank = 0){
    if(all(class(input) == "character")){
        coord <- splitCoords(input)
        coord$start <- coord$start - flank
        coord$end <- coord$end + flank
        coord <- paste0(coord$chr, ":", coord$start, "-", coord$end)    
        return(coord)
    }
    if( "data.frame" %in% class(input)){
        stopifnot( all(c("chr", "start", "end") %in% names(input) ) | all(c("chr","pos") %in% names(input) )  )
        stopifnot( flank >= 0)
        if( all( c("chr","start","end") %in% names(input) ) ){
            coords <- paste0(input$chr, ":", input$start - flank, "-",input$end + flank) 
        }
        if( all( c("chr", "pos") %in% names(input) ) ){
            coords <- paste0(input$chr, ":", input$pos - (flank + 1), "-", input$pos + flank)
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

extractLoci <- function(gwas){
    loci_path <- gwas$top_filtered_path
    
    # stupid exception for Daner_2020 - filtered loci are totally different file right now
    if(gwas$dataset == "Daner_2020"){
        loci_path <- gwas$top_path
    }
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
    
    stopifnot( all( c(gwas$top_locus, gwas$top_snp, gwas$top_chrom) %in% colnames(loci_df) ) )
    loci_df$locus <- loci_df[[gwas$top_locus]] 
    loci_df$snp <- loci_df[[gwas$top_snp]]
    loci_df$chr <- loci_df[[gwas$top_chrom]]
    loci_df$pos <- as.numeric(loci_df[[gwas$top_pos]])
    
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
    stopifnot( "chr" %in% names(data) )
    #stopifnot( all(!is.na(data$snp) ))
    chr <- gsub("chr", "", unique(data$chr))
    stopifnot( !is.na(chr) )
    stopifnot( length(chr) == 1)
    # subset maf by chromosome
    maf <- maf[chr] 
    # match on RSID
    message("   * before joining MAF: ", nrow(data), " SNPs")
    data$MAF <- maf$MAF[match(data$snp, maf$snp)]
    matched <- data[ !is.na(data$MAF) ,]
    message("   * after joining MAF: ", nrow(matched), " SNPs")
    if(nrow(matched) == 0 ){
        print(head(data) )
    stop("no SNPs match on RSID. Inspect data")
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
    gwas_path <- file.path( refFolder,  paste0(gwas$dataset, ".processed.tsv.gz" ))
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
     
    cmd <- paste( "ml bcftools; tabix -h ", gwas_path, coord ) 
    message(" * running command: ", cmd)
 
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    message(" * GWAS extracted!")
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
        message("       * MAF not present - using 1000 Genomes MAF")
        # how to get this in to the function? it can't find maf_1000g
        result <- matchMAF(result, maf = maf_1000g)
    }else{
        message("       * using supplied MAF")
        names(result)[names(col_dict) == mafCol] <- "MAF" 
    }

    # convert standard error to the variance
    # se is standard error - not standard deviation!
    result$varbeta <- result$varbeta^2
    
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
    message(" * lifting over....")
    lifted_over <- rtracklayer::liftOver(coord_gr, chain )
    message(" * lifted!" )
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
    
    message( paste0( nrow(sig), " significant genes or splicing events at this locus" ) )
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
    cmd <- paste( "ml bcftools; tabix -h ", qtl$full_path, coord )
    message("       * running command: ", cmd)
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    # deal with regions of no QTL association
    if( nrow(result) == 0){ 
        return(NULL) 
    }
    #stopifnot( nrow(result) > 0 )
    message("       * QTL extracted!")
    return(result)
}


# extract all nominal QTL P-values overlapping the flanked GWAS hit
# split by Gene being tested
# sig_level doesn't do anything yet
extractQTL <- function(qtl, coord, sig_level = 0.05, force_maf = FALSE){
    # variables stored in qtl
    # if qtl type is parquet then read in parquet files
    # if qtl type is tabix then query the coordinates through tabix 
    

    # check columns
    stopifnot( all(!is.na(c(qtl$full_snp, qtl$full_pheno, qtl$full_p, qtl$full_effect, qtl$full_se, qtl$N, qtl$build) ) ))
    pvalCol <- qtl$full_p
    betaCol <- qtl$full_effect
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
    stopifnot( ncol(result) == length(col_dict) )
    names(result)[names(col_dict) == phenoCol] <- "gene"
    names(result)[names(col_dict) == mafCol]   <- "MAF"
    names(result)[names(col_dict) == pvalCol]  <- "pvalues"
    names(result)[names(col_dict) == betaCol]  <- "beta"
    names(result)[names(col_dict) == seCol]    <- "varbeta"
    names(result)[names(col_dict) == snpCol]   <- "snp"
    # Young microglia use log10P - convert
    if( pvalCol == "log10_p"){
        result$pvalues <- 10^result$pvalues
    }

    # deal with MAF - meta-analysis outputs won't have it
    # this requires "chr" to be present in QTL data
    if( is.na(mafCol) | force_maf == TRUE ){
        message("       * MAF not present - using 1000 Genomes MAF")
        #stopifnot( "chr" %in% names(col_dict) )
        names(result)[names(col_dict) == "chr"] <- "chr"
        result <- matchMAF(result, maf = maf_1000g)
    }else{
        message("       * using supplied MAF")
        names(result)[names(col_dict) == mafCol] <- "MAF"
    }

    # deal with edge cases - 1 gene and the gene id is NA
    if( all( is.na(result$gene) ) ){
        message(" * no genes in QTL result")
        return(NULL)
    }else{
        result <- dplyr::filter(result, !is.na(gene) )
    }
    # if MAF needs matching then chrCol will be "chr" - change it back to "QTL_chr"
    names(result)[names(col_dict) == chrCol]   <- "QTL_chr"
    names(result)[names(col_dict) == posCol]   <- "QTL_pos"

    #print(head(result) ) 
    
    # don't forget to square the standard error to get the variance             
    result$varbeta <- result$varbeta^2
 
    # remove rows with variance of 0 - this will corrupt COLOC
    result <- dplyr::filter(result, varbeta != 0) 

    # retain only associations within locus coords
    res_subset <- dplyr::select(result, gene, snp, pvalues, beta, varbeta, MAF, QTL_chr, QTL_pos) #%>%
        #dplyr::filter( pos >= coord_split$start & pos <= coord_split$end)

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

       


# running COLOC between a GWAS locus and all eQTLs within 1MB either side
# flank coordinates
# extract GWAS nominal SNPs from region
# find GWAS hit SNP for recording
# extract QTLs within range - at set significance level
# join together and run COLOC
# output a list of objects:
## COLOC results - this should include the GWAS locus, the gene of interest, the top QTL variant and the top QTL p-value
## object - this should be a table combining the two input datasets, inner joined on "snp", to be used for plotting.
runCOLOC <- function(gwas, qtl, hit){
    message(" * Analysing LOCUS: ", hit$locus )
    # hit is a dataframe containing "snp", "chr", "pos", "locus"
    
    hit_coord <- joinCoords(hit, flank = 0)
    hit_range <- joinCoords(hit, flank = 1e6)
 
    # extract GWAS and QTL summary for given coord
    message("       * extracting GWAS")
    g <- extractGWAS(gwas, hit_range)
    # get hit info from GWAS
    hit_snp <- hit$snp
    # this assumes that the hit SNP is within the g object
    # if the top hit doesn't have an RSID or cannot be MAF matched then this is violated
    # I will have to use the hit object instead and accept that not all top_loci lists include the MAF, P value etc.
    #if( hit_snp == hit$snp ){
    #hit_info <- as.data.frame(g, stringsAsFactors = FALSE) %>% filter(snp == hit_snp) %>% 
    #            select(GWAS_SNP = snp, GWAS_P = pvalues, GWAS_effect_size = beta, GWAS_MAF = MAF) %>%
    #            mutate(locus = hit$locus) %>% select(locus, everything() )
    #}else{
    # now just use top loci
   
    # problem when no p column in hit!

    hit_info <- hit %>%
            select(locus, GWAS_SNP = snp, GWAS_P = p)
    #}
 
    # if GWAS is hg19 then lift over to hg38
    # use the lifted over position in hit_info
    if(gwas$build != qtl$build){
        qtl_coord <- liftOverCoord(hit_coord, from = gwas$build, to = qtl$build)
    }else{
        qtl_coord <- hit_coord
    }
    # record the position of the GWAS top SNP
    hit_info$GWAS_chr <- splitCoords(qtl_coord)$chr
    hit_info$GWAS_pos <- splitCoords(qtl_coord)$end


    qtl_range <- joinCoords(qtl_coord, flank = 1e6)
    
    message("       * extracting QTLs")
    q <- extractQTL(qtl, qtl_range)
    if( is.null(q) ){ return(NULL) }
    # for each gene extract top QTL SNP    
    qtl_info <- q %>% 
        map_df( ~{ 
            as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_Beta = beta, QTL_MAF = MAF, QTL_chr, QTL_pos)
        }) 
    
    # actually run COLOC
        message("      * running COLOC")
    # run COLOC on each QTL gene 
    # return the results table and an object containing g and q
    coloc_res <- 
        purrr::map( q, ~{
            # if QTL has no overlapping SNPs with GWAS summary then return NULL
            if( length(intersect(g$snp, .x$snp) ) == 0 ){ 
                warning("GWAS and QTL have no SNPs in common")
                return(NULL) 
            }
            coloc_object <- coloc.abf(dataset1 = g, dataset2 = .x)
            # add in g and q to coloc object
            gq <- dplyr::inner_join(as.data.frame(g),as.data.frame(.x), by = "snp", suffix = c(".gwas", ".qtl") )
            coloc_object$results <- dplyr::inner_join(coloc_object$results, gq, by = "snp" )
            # pull out coloc summary as data.frame
            coloc_df <- as.data.frame(t(coloc_object$summary), stringsAsFactors = FALSE)
            return( list(df = coloc_df, object = coloc_object) )
        })
    message("       * COLOC finished")
    coloc_df <- map_df(coloc_res, "df", .id = "gene") %>% 
        dplyr::mutate(locus = hit$locus ) %>%
        dplyr::select(locus, gene, everything() )
    message("       * COLOC res created")
    # bind coloc_res to other tables
    full_res <- full_join(qtl_info, coloc_df, by = "gene")
    full_res <- full_join(hit_info, full_res, by = "locus")

    # add pairwise LD between GWAS and QTL SNP
    #full_res <- calc_LD(full_res)

    return(list(full_res = full_res, coloc_object = coloc_res) ) 
}

# for each COLOC result per locus
# calculate pairwise LD between GWAS SNP and QTL top SNPs
calc_LD <- function( coloc_res ){
    require(LDlinkR)
    # get LDlink token from .Renviron (only Jack has this, for access go to https://ldlink.nci.nih.gov/?tab=apiaccess)
    token <- Sys.getenv("LDLINK_TOKEN")
    if( token == ""){ 
        warning(" * LDlink token not found!" )
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
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = ""),
        make_option(c('--gwas', '-g'), help= "the dataset ID for a GWAS in the GWAS/QTL database" ),
        make_option(c('--qtl', '-q'), help = "the dataset ID for a QTL dataset in the GWAS/QTL database"),
        make_option(c('--debug'), help = "load all files and the nsave RData without running COLOC", action = "store_true", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFolder <- opt$outFolder
gwas_dataset <- opt$gwas
qtl_dataset <- opt$qtl
debug <- opt$debug
#gwas_prefix <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/processed/Nicolas_2018_processed_"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_expression/peer30/LumbarSpinalCord_expression_peer30"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_splicing/peer20/LumbarSpinalCord_splicing_peer20"

#hits_file <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/COLOC/Nicolas_2018/Nicolas_2018_hits_1e-7.tsv"
#outFile <- "test_coloc_results.tsv"

#options(echo = TRUE)

# load in data
# extract top GWAS loci
# for each locus run COLOC

library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

maf_1000gp1 <- "/sc/arion/projects/ad-omics/data/references/1KGP1/1000G_EUR_MAF.bed.gz"
maf_1000gp3 <- "/sc/arion/projects/ad-omics/data/references/1KGPp3v5/EUR_MAF/EUR.all.phase3_MAF.bed.gz"

# load in MAF table

if( !exists("maf_1000g")){

maf_1000g <- loadMAF(maf_1000gp3)
# load in liftover chain
chain_hg19_hg38 <- import.chain("/sc/arion/projects/ad-omics/data/references/liftOver/hg19ToHg38.over.chain")    

}

main <- function(){
    
        #gwas_dataset <- "Kunkle_2019"
    #gwas_dataset <- "Nalls23andMe_2019"
    #gwas_dataset <- "Marioni_2018"
    #qtl_dataset <- "Microglia_THA"
    
    gwas <- pullData(gwas_dataset, "GWAS")
    qtl <- pullData(qtl_dataset, "QTL")
    
    top_loci <- extractLoci(gwas)
    
    # for testing
    #top_loci <- top_loci[2,]
    if( debug == TRUE){ save.image("debug.RData") }
  
    all_coloc <- purrr::map(
        1:nrow(top_loci), ~{
            res <- runCOLOC(gwas, qtl, hit = top_loci[.x,])
            if(is.null(res) ){return(NULL) }
            return(res)
        }
    )
    all_res <- map_df(all_coloc, "full_res")
    all_obj <- map(all_coloc, "coloc_object")    
    names(all_obj) <- top_loci$locus
    
    
    # arrange by H4
    all_res <- arrange(all_res, desc(PP.H4.abf)  )
    
    outFile <- paste0(outFolder, qtl_dataset, "_", gwas_dataset, "_COLOC.tsv")
    readr::write_tsv(all_res, path = outFile)
    # save COLOC objects for plotting
    save(all_obj, file = gsub("tsv", "RData", outFile))
}

main()


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
# run MR to get different Medelian Randomisation.
# extract top QTL SNP for locus
# record results in table - GWAS locus, Gene, top QTL SNP, MR p value


suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(optparse))

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

option_list <- list(
  make_option(c('-o', '--outFolder'), help='the path to the output file', default = ""),
  make_option(c('--gwas', '-g'), help= "the dataset ID for a GWAS in the GWAS/QTL database" ),
  make_option(c('--qtl', '-q'), help = "the dataset ID for a QTL dataset in the GWAS/QTL database"),
  make_option(c('--fdr', '-f'), help = "the exposure filter option, if TRUE, use FDR < 0.05", default = FALSE),
  make_option(c('--cores', '-c'), help = "the number of cores for parallelizing", default = 1),
  make_option(c('--plink_bin'), help = "Path to plink binary (optional)", default = ""),
  make_option(c('--debug'), help = "load all files and then save RData without running MR", action = "store_true", default = FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFolder <- opt$outFolder
gwas_dataset <- opt$gwas
qtl_dataset <- opt$qtl
fdr_filter <- opt$fdr
debug <- opt$debug
ncores <- opt$cores
plink_bin <- opt$plink_bin

if (!is.null(plink_bin) && nzchar(plink_bin) && file.exists(plink_bin)) {
  message("Using plink from config: ", plink_bin)
  PLINK_BIN <- plink_bin
} else {
  message("Using genetics.binaRies::get_plink_binary()")
  PLINK_BIN <- genetics.binaRies::get_plink_binary()
}

# flank coordinates by set number of bases (default 1MB)
# work on either a coordinate string or a dataframe containing chr start and end columns
joinCoords <- function(input, flank = 0){
  if(all(class(input) == "character")){
    coord <- splitCoords(input)
    coord$start <- coord$start - flank
    coord$end <- coord$end + flank
    if(coord$start < 0 ){
      coord <- paste0(coord$chr, ":", "0", "-", coord$end)    
    }else{
      coord <- paste0(coord$chr, ":", coord$start, "-", coord$end)    
    }
    return(coord)
  }
  if( "data.frame" %in% class(input)){
    stopifnot( all(c("chr", "start", "end") %in% names(input) ) | all(c("chr","pos") %in% names(input) )  )
    stopifnot( flank >= 0 )
    if( all( c("chr","start","end") %in% names(input) ) ){
      if((as.numeric(input$start) - flank) < 0 ){
        coords <- paste0(input$chr, ":", "0", "-",input$end + flank) 
      }else{
        coords <- paste0(input$chr, ":", input$start - flank, "-",input$end + flank) 
      }
    }
    if( all( c("chr", "pos") %in% names(input) ) ){
      if((as.numeric(input$pos) - (flank + 1)) < 0 ){
        coords <- paste0(input$chr, ":", "0", "-", input$pos + flank)
      }else{
        coords <- paste0(input$chr, ":", input$pos - (flank + 1), "-", input$pos + flank) 
      }
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
    message(" * loading in MAF from gnomAD v4")
    maf <- data.table::fread(path, nThread = 4)
    names(maf) <- c("chr", "pos", "snp", "MAF")
    maf$chr <- gsub('chr', '', maf$chr)
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
     
    cmd <- paste( "ml tabix; tabix -h ", gwas_path, coord ) 
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
    if( debug == TRUE){ result$MAF <- NA }else{
        if( is.na(mafCol) | force_maf == TRUE ){
            message("       * MAF not present - using 1000 Genomes MAF")
            # how to get this in to the function? it can't find maf_1000g
            result <- matchMAF(result, maf = maf_1000g)
        }else{
            message("       * using supplied MAF")
            names(result)[names(col_dict) == mafCol] <- "MAF" 
        }
    }
    # convert standard error to the variance
    # se is standard error - not standard deviation!
    result$varbeta <- result$varbeta^2
    print(head(result))
    #save.image("debug.RData") 
    result <- dplyr::select( result, snp, pvalues, beta, varbeta, MAF, chr, pos, A1, A2)
    
    #result <- as.list(result)
    
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

extractQTL_top <- function(qtl, fdr_filter, sig_level = 0.05){

  if(!file.exists(qtl$top_path)){ return(NULL) }
  
  # cmd <- paste( "ml bcftools; tabix -h ", qtl$full_path, coord )
  message("       * Loading the Top QTL: ", qtl$top_path)
  result <- as.data.frame(data.table::fread(qtl$top_path, nThread = 4) )
  # deal with regions of no QTL association
  if(is.null(result) ){ return(NULL) }
  
  # check columns
  stopifnot( all(!is.na(c(qtl$top_snp, qtl$top_pheno, qtl$top_p, qtl$N, qtl$build) ) ))
  pvalCol <- qtl$top_p
  betaCol <- qtl$top_effect # if not present will be set to NA
  phenoCol <- qtl$top_pheno
  seCol <- qtl$top_se
  snpCol <- qtl$top_snp
  mafCol <- qtl$top_maf
  chrCol <- qtl$top_chrom
  posCol <- qtl$top_pos
  a1Col <- qtl$top_A1
  a2Col <- qtl$top_A2
  ppermCol <- qtl$top_p_perm
  qCol <- qtl$top_q
  
  
  # assign column names
  cmd <- paste( "zless ", qtl$top_path, " | head -1 " )
  columns <- colnames(data.table::fread( cmd = cmd ))
  col_dict <- setNames(1:length(columns), columns)
  
  print(col_dict)
  #stopifnot( ncol(result) == length(col_dict) )
  col_dict <- col_dict[1:ncol(result)] # hacky workaround for now
  names(result)[names(col_dict) == phenoCol] <- "gene"
  names(result)[names(col_dict) == mafCol]   <- "MAF"
  names(result)[names(col_dict) == pvalCol]  <- "pvalues"
  names(result)[names(col_dict) == a1Col] <- "A1"
  names(result)[names(col_dict) == a2Col] <- "A2"
  names(result)[names(col_dict) == ppermCol] <- "pval_perm"
  names(result)[names(col_dict) == qCol] <- "qvalue"
  
  if("qvalue" %in% colnames(result)){
    result <- dplyr::filter(result, qvalue < sig_level) 
    message(Sys.time()," * QTL qvalue filter : ", nrow(result))
  }
  # if(fdr_filter){
  #   result <- dplyr::filter(result, pval_perm < 0.05)
  #     message(Sys.time()," * QTL FDR filter : ", nrow(result))
  #   }else{
  #     result <- dplyr::filter(result, pvalues < 5e-08)
  #   message(Sys.time()," * QTL Filter 5e-08 filter : ", nrow(result))
  # }
  # 
  if(nrow(result)==0){
    message(Sys.time()," * No significant QTL")
    return(NULL)
    }
  
  # print(head(result))
  genelist <- result$gene
  message("       * Sig QTLs Gene count : ", length(genelist))
  #stopifnot( nrow(result) > 0 )
  
  message("       * QTL top extracted!")
  
  return(genelist)
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
    cmd <- paste( "ml tabix; tabix -h ", qtl$full_path, coord )
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
    stopifnot( all(!is.na(c(qtl$full_snp, qtl$full_pheno, qtl$full_p, qtl$N, qtl$build) ) ))
    pvalCol <- qtl$full_p
    betaCol <- qtl$full_effect # if not present will be set to NA
    phenoCol <- qtl$full_pheno
    seCol <- qtl$full_se
    snpCol <- qtl$full_snp
    mafCol <- qtl$full_maf
    chrCol <- qtl$full_chrom
    posCol <- qtl$full_pos
    a1Col <- qtl$full_A1
    a2Col <- qtl$full_A2
    fdrCol <- qtl$full_p_perm
    
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
    
    print(col_dict)
    #stopifnot( ncol(result) == length(col_dict) )
    col_dict <- col_dict[1:ncol(result)] # hacky workaround for now
    names(result)[names(col_dict) == phenoCol] <- "gene"
    names(result)[names(col_dict) == mafCol]   <- "MAF"
    names(result)[names(col_dict) == pvalCol]  <- "pvalues"
    names(result)[names(col_dict) == a1Col] <- "A1"
    names(result)[names(col_dict) == a2Col] <- "A2"
    names(result)[names(col_dict) == fdrCol] <- "FDR"
    
    if( !is.na(betaCol) ){
      names(result)[names(col_dict) == betaCol]  <- "beta"
    }
    if( !is.na(seCol) ){
      names(result)[names(col_dict) == seCol]    <- "varbeta"
      # don't forget to square the standard error to get the variance             
      result$varbeta <- result$varbeta^2
      
      # remove rows with variance of 0 - this will corrupt COLOC
      result <- dplyr::filter(result, varbeta != 0) 
    }
    names(result)[names(col_dict) == snpCol]   <- "snp"
    # Young microglia use log10P - convert
    if( pvalCol == "log10_p"){
      result$pvalues <- 10^result$pvalues
    }
    
    # deal with MAF - meta-analysis outputs won't have it
    # this requires "chr" to be present in QTL data
    if( debug == TRUE){ result$MAF <- NA }else{
      if( is.na(mafCol) | force_maf == TRUE ){
        message("       * MAF not present - using 1000 Genomes MAF")
        #stopifnot( "chr" %in% names(col_dict) )
        names(result)[names(col_dict) == "chr"] <- "chr"
        result <- matchMAF(result, maf = maf_1000g)
        print(head(result) )
      }else{
        message("       * using supplied MAF")
        names(result)[names(col_dict) == mafCol] <- "MAF"
        print(col_dict)
        print(head(result) )
      }
    }
    # deal with edge cases - 1 gene and the gene id is NA
    if( all( is.na(result$gene) ) ){
      message(" * no genes in QTL result")
      return(NULL)
    }else{
      result <- dplyr::filter(result, !is.na(gene) )
    }
    # if MAF needs matching then chrCol will be "chr" - change it back to "QTL_chr"
    print( names(col_dict) )
    print( names(col_dict) == chrCol )
    names(result)[which(names(col_dict) == chrCol) ]   <- "QTL_chr"
    names(result)[which(names(col_dict) == posCol) ]   <- "QTL_pos"
    
    #print(head(result) ) 
    
    # retain only associations within locus coords
    if( !is.na(betaCol) & !is.na(seCol) &  !is.na(fdrCol)){
      res_subset <- dplyr::select(result, gene, snp, pvalues, beta, varbeta, MAF, QTL_chr, QTL_pos, A1, A2, FDR)
    }else{
      res_subset <- dplyr::select(result, gene, snp, pvalues, MAF, QTL_chr, QTL_pos, A1, A2)
    }
    print(head(res_subset) ) 
    print("passed this point")
       
    # split by gene, convert to list object

    return(res_subset)
}

       
MR_format <- function(dataset , types){
  format_data(
    dataset,
    type = types,
    snps = NULL,
    header = TRUE,
    phenotype_col = types,
    effect_allele_col = "A1",
    other_allele_col = "A2",
    snp_col = "snp",
    beta_col = "beta",
    se_col = "varbeta",
    eaf_col = "MAF",
    pval_col = "pvalues"
  ) ##eaf_col not necessary here?
}

# running MR between a GWAS locus and all eQTLs within 1MB either side
# flank coordinates
# extract GWAS nominal SNPs from region
# find GWAS hit SNP for recording
# extract QTLs within range - at set significance level
# join together and run MR
# output a list of objects:
## MR results - this should include the GWAS locus, the gene of interest, the top QTL variant and the top QTL p-value
## object - this should be a table combining the two input datasets, inner joined on "snp", to be used for plotting.
runMR <- function(gwas, qtl,  hit, sig_gene, fdr_filter){

    message(" * Analysing LOCUS: ", hit$locus )
    # hit is a dataframe containing "snp", "chr", "pos", "locus"
    
    hit_coord <- joinCoords(hit, flank = 0)
    hit_range <- joinCoords(hit, flank = 1e6)
 
    # extract GWAS and QTL summary for given coord
    message("       * extracting GWAS")
    g <- extractGWAS(gwas, hit_range)
    g$outcome <- c(hit$locus)
    message(Sys.time()," * GWAS raw : ", nrow(g))
    
    gwas_snplist <- g$snp

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
    # if(!is.null(sig_gene)){
    #   sig_gene <- sapply(strsplit(sig_gene, '\\.'), '[[', 1)
    # }

    # coloc_info <- coloc %>% dplyr::filter((coloc$locus==hit$locus))
    
    if(is.null(q)){
      message("       * No significant gene!")
      return(NULL)
    }
  
    
    q <- q[which(q$snp %in% gwas_snplist),]
    if( nrow(q)==0){ 
      message("       * No overlap eQTL with GWAS")
      return(NULL) 
    }   
    message(Sys.time()," * Overlap eQTL with GWAS : ", nrow(q))
    # split by gene, convert to list object
    q_info <- 
      split(q, q$gene) %>%
      purrr::map( ~{ 
        x = as.list(.x)
        return(x)
      })
    # for each gene extract top QTL SNP
    if( all( c("beta", "varbeta","FDR") %in% names(q[[1]]) ) ){
      qtl_info <- q_info %>% 
        map_df( ~{ 
          as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_Beta = beta, QTL_SE = varbeta, QTL_MAF = MAF, QTL_chr, QTL_pos) %>%
            mutate( QTL_SE = sqrt(QTL_SE) ) # get SE back from Variance
        })
      # get out beta, se and P for GWAS SNP 
      gwas_info <- q_info %>%
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
      print(head(qtl_info) )
    }else{
      qtl_info <- q_info %>% 
        map_df( ~{ 
          as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_MAF = MAF, QTL_chr, QTL_pos)
        }) 
      
    }
    rm(q_info)


    # actually run MR
    message("      * running MR")
    
    g_mr <- MR_format(g, 'outcome')
    
    MR = list()
    genelist <- q$gene %>% unique()
    message("       * QTLs Gene count : ", length(genelist))
    genelist <- genelist[which(genelist %in% sig_gene)]
    message("       * QTLs Gene count : ", length(genelist))
    if( length(genelist)==0){ 
      message("       * No sig QTLs Gene")
      return(NULL) 
    }   
    message("       * QTLs Gene list : ", genelist)
  
    for(i in 1:length(genelist)){
         genename <- genelist[i]
         message("       * QTLs Gene : ", genename)
          
         # coloc_info <- coloc %>% dplyr::filter((coloc$QTL_Ensembl==genename)&(coloc$locus==hit$locus))
         
         message(Sys.time()," * QTL raw : ", nrow(q))
         q_tmp <- q %>% dplyr::filter(gene == genename)
         
         
         message(Sys.time()," * QTL ",genename," : ", nrow(q_tmp))
         if(nrow(q_tmp)==0){next}
         q_tmp$exposure <- genename
         
         if(fdr_filter){
           if(("FDR" %in% colnames(q_tmp))&(genename %in% sig_gene)){
             q_tmp <- dplyr::filter(q_tmp, FDR < 0.05)
             message(Sys.time()," * QTL FDR filter : ", nrow(q_tmp))
           }else{
             message("       * No significant eQTL")
             next
           }
         }else{
           q_tmp <- dplyr::filter(q_tmp, pvalues < 5e-08)
           message(Sys.time()," * QTL Filter 5e-08 filter : ", nrow(q_tmp))
         }
    
         
         if(nrow(q_tmp)==0){
           message("       * No significant eQTL")
           next
           }
         
         # q_tmp <- arrange(q_tmp, FDR )
         # print(head(q_tmp))
         
         q_tmp <- arrange(q_tmp, pvalues )
         print(head(q_tmp))
         
        names(q_tmp)[which(names(q_tmp) == 'snp') ]   <- "rsid"
        names(q_tmp)[which(names(q_tmp) == 'pvalues') ]   <- "pval" 
        
        
        message("       * Clumping With LD clump")
        try_ld_clump <- try(
          #(S_g <- runsusie(g, n = n_g, maxit = 100, repeat_until_convergence = F) ), silent=FALSE)
          (q_clump <- ieugwasr::ld_clump(q_tmp,
                                         clump_kb = 1000,
                                         clump_r2 =  0.01,
                                         # pop = 'EUR', 
                                         # opengwas_jwt = Sys.getenv('LDCLUMP_TOKEN'),
                                         plink_bin = PLINK_BIN,                                         
                                         bfile ='/sc/arion/projects/ad-omics/data/references/LDreference/EUR'
           ) ), silent=FALSE)
         if (class(try_ld_clump) == "try-error"){
           # This error happens when ld clump was not converged
           message("       * No Clumping With LD clump")
           next
         }
         print(head(q_clump))
        
    
      
         if(nrow(q_clump)==0){next}
         q_clump$pval <- as.numeric(q_clump$pval)
         q_clump$beta <- as.numeric(q_clump$beta)
         q_clump$varbeta <- as.numeric(q_clump$varbeta)
        
         names(q_clump)[which(names(q_clump) == 'rsid') ]   <- "snp"
         names(q_clump)[which(names(q_clump) == 'pval') ]   <- "pvalues"
         
         
         message("       * Formatting")
         q_mr <- MR_format(q_clump,'exposure')
         
         message("       * Harmonisation")
        
         dat <- harmonise_data(
             exposure_dat = q_mr,
             outcome_dat = g_mr 
         )
         if (nrow(dat) == 0){
           message("       * None of the provided SNPs can be used for MR analysis")
           next
         }
         
         message("       * Running MR")
         
         try_mr <- try(
           #(S_g <- runsusie(g, n = n_g, maxit = 100, repeat_until_convergence = F) ), silent=FALSE)
           (res <- mr(dat) ), silent=FALSE
           )
         if (class(try_mr) == "try-error"){
           # This error happens when ld clump was not converged
           message("       * No MR result")
           next
         }
         if(is.null(res$b)){next} 
         # print(head(res))
    
         generate_odds_ratios(res) 
    #      het <-mr_heterogeneity(dat)
    #      plt <- mr_pleiotropy_test(dat)
    #      sin <- mr_singlesnp(dat)
    #      result <- combine_all_mrresults(res,
    #                                      het,
    #                                      plt,
    #                                      sin
    # )
    
         res$locus <- res$outcome
         print(res)
         # mr_results <- dplyr::inner_join(as.data.frame(coloc_info),as.data.frame(res), by = "locus" ) 
         MR[[i]]=res
        
    }
    message("       * MR finisehd")
    if( all(sapply(MR, is.null))){ return(NULL) }
    full_res <- as.data.frame(do.call(rbind,MR))
    
    print(head(full_res))
    full_res$gene <- full_res$exposure
    full_res <- right_join(qtl_info, full_res, by = "gene")
    full_res <- full_join(hit_info, full_res, by = "locus")
    
    return(full_res) 
}

# load in data
# extract top GWAS loci
# for each locus run COLOC
maf_1000gp1 <- "/sc/arion/projects/ad-omics/data/references/1KGP1/1000G_EUR_MAF.bed.gz"
maf_1000gp3 <- "/sc/arion/projects/ad-omics/data/references/1KGPp3v5/EUR_MAF/EUR.all.phase3_MAF.bed.gz"
gnomad_v4 <- '/sc/arion/projects/ad-omics/data/references/gnomAD_v4/genomes/EUR_MAF/EUR.all.gnomAD_v4_MAF0001.bed.gz'

# load in MAF table
if( !exists("maf_1000g")){

        maf_1000g <- loadMAF(gnomad_v4)
        # load in liftover chain
        chain_hg19_hg38 <- import.chain("/sc/arion/projects/ad-omics/data/references/liftOver/hg19ToHg38.over.chain")    
    }

main <- function(){
    
    if(fdr_filter){ 
      message("       * MR filter option : FDR")
    }else{
      message("       * MR filter option : GWAS p-value")
    }
  
    gwas <- pullData(gwas_dataset, "GWAS")
    qtl <- pullData(qtl_dataset, "QTL")
    
    sig_gene <- extractQTL_top(qtl, fdr_filter)
    # stopifnot(is.null(sig_gene))

      top_loci <- extractLoci(gwas)
      
      if( debug == TRUE){ save.image("debug.RData") }
      
      MR_results =list()
      
      MR_results <- mclapply(seq_len(nrow(top_loci)), function(i){
         runMR(gwas, qtl, hit = top_loci[i,], sig_gene, fdr_filter)}, mc.cores = ncores)

      if( all(sapply(MR_results, is.null))){
          if(fdr_filter){
            all_res <- data.frame()
            outFile <- paste0(outFolder, qtl_dataset, "_", gwas_dataset, "_MR_wFDR.tsv")
            message("       * MR finished : ", outFile)
            readr::write_tsv(all_res, file = outFile)
    
          }else{
            all_res <- data.frame()
            outFile <- paste0(outFolder, qtl_dataset, "_", gwas_dataset, "_MR.tsv")
            message("       * MR finished : ", outFile)
            readr::write_tsv(all_res, file = outFile)
          }
      }else{
        MR_results = as.data.frame(do.call(rbind, MR_results))
        
        MR_results <- arrange(MR_results, desc(pval)  )
        print(head(MR_results))
        if(fdr_filter){
          outFile <- paste0(outFolder, qtl_dataset, "_", gwas_dataset, "_MR_wFDR.tsv")
          message("       * MR finished : ", outFile)
          readr::write_tsv(MR_results, file = outFile)
          
        }else{
          outFile <- paste0(outFolder, qtl_dataset, "_", gwas_dataset, "_MR.tsv")
          message("       * MR finished : ", outFile)
          readr::write_tsv(MR_results, file = outFile)
        }
      }
}
    
main()


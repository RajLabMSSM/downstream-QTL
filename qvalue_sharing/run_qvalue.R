
# Q-value sharing using the GWAS/QTL database
# requires development version of qvalue
# remotes::install_github("jdstorey/qvalue")

# get QTL metadata for source (top QTLs) and target (all QTLs) from database
pullData <- function(dataset, type = "GWAS"){
    message(Sys.time()," * selected dataset: ", dataset)
    db_path <- "/sc/hydra/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"

    message(Sys.time()," * reading GWAS database from ", db_path)
    stopifnot( file.exists(db_path) )
    if( type == "GWAS"){
        n_sheet <- 2
    }
    if( type == "QTL"){
        n_sheet <- 3
    }
    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = n_sheet,na= c("", "-","NA")))
    
    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )

    #stopifnot( all( c("full_path", "full_chrom", "full_pos") %in% names(gwas) ) )

    #stopifnot( file.exists(gwas$full_path) )

    return(gwas)
}
# read in top permutation QTL results
extractTopQTL <- function(qtl, qvalue_threshold){
    loci_path <- qtl$top_path
    # loci are stored either in comma-separated (csv), tab-separated files (txt, tsv) or in excel files (xlsx, xlsm)
    # read in loci file depending on file type
    message(" * reading in top QTLs from ", loci_path)
    loci_df <- readr::read_tsv(loci_path)
    
    stopifnot( nrow(loci_df) > 1) 
    # rename columns
    # return only core set of columns - locus, SNP, q-value, phenotype (Gene or splicing event)
    stopifnot( all( c("top_q", "top_p_perm", "top_snp", "top_pheno", "top_snp") %in% names(qtl) ) ) 
    stopifnot( all(!is.na(c(qtl$top_snp, qtl$top_pheno, qtl$top_p_perm) ) ) )
    loci_df$snp <- loci_df[[qtl$top_snp]]
    loci_df$pval <- loci_df[[qtl$top_p_perm]]
    loci_df$pheno <- loci_df[[qtl$top_pheno]]
    
    # adjust EnsemblID - remove tag number
    loci_df$pheno <- stringr::str_split_fixed(loci_df$pheno, "\\.", 2)[,1]
    
    # calculate q-value if not available
    if( !is.na(qtl$top_q) ){
        loci_df$qval <- loci_df[[qtl$top_q]]
    }else{
        loci_df$qval <- qvalue::qvalue(loci_df$pval)$qvalues 
    }
    
    # use only significant associations at given q-value
    loci_df <- dplyr::filter(loci_df, qval < qvalue_threshold )
     
    # if sQTL then process the sQTL phenotype     
 
    loci_clean <- dplyr::select(loci_df, snp, pheno, pval, qval)
    # same for beta and MAF?
    return(loci_clean)
}

# for top QTLs, read in per chromosome and extract positions
extractTargetQTL <- function(qtl, chr, source_qtls){
    require(dplyr)
    stopifnot( file.exists(qtl$full_path) )
     
    # read in QTL
    cmd <- paste0("ml bcftools; tabix ", qtl$full_path, " ", chr)
    message( " * ", cmd) 
    result <- data.table::fread( cmd = cmd, nThread = 8 )
    if( nrow(result) == 0){ 
        return(NULL) 
    }

    # rename columns 
     stopifnot( all(!is.na(c(qtl$full_snp, qtl$full_pheno, qtl$full_p, qtl$full_beta, qtl$full_se, qtl$full_maf, qtl$N, qtl$build) ) ))
    pvalCol <- qtl$full_p
    betaCol <- qtl$full_beta
    phenoCol <- qtl$full_pheno
    seCol <- qtl$full_se
    snpCol <- qtl$full_snp
    mafCol <- qtl$full_maf
    if(is.null(result) ){ return(NULL) }

    cmd <- paste( "zless ", qtl$full_path, " | head -1 " )
    columns <- colnames(data.table::fread( cmd = cmd, header= TRUE))
    col_dict <- setNames(1:length(columns), columns)
    #stopifnot( ncol(result) == length(col_dict) )
    names(result)[names(col_dict) == phenoCol] <- "pheno"
    names(result)[names(col_dict) == pvalCol]  <- "pvalue"
    names(result)[names(col_dict) == snpCol]   <- "snp"
    
    # weird exception - Young Microglia stores p-values as log10 - convert
    if( qtl$full_p == "log10_p"){
        result$pvalue <- 10^result$pvalue
    }
    #return(result)
   
    # adjust EnsemblID - remove tag number
    result$pheno <- stringr::str_split_fixed(result$pheno, "\\.", 2)[,1]

 
    # subset out just SNP-Gene pairs in source_qtls
    res_subset <- dplyr::filter(result, pheno %in% source_qtls$pheno & snp %in% source_qtls$snp ) %>%
        select( pheno, snp, pvalue_target = pvalue)
    
    return(res_subset)
}

# run qvalue on set of P-values
# handle errors
run_qvalue <- function(x){
    tryCatch(
        expr = {
            return(qvalue::qvalue_truncp(x) )
            message("Successfully executed the log(x) call.")
        },
        error = function(e){
            message('Caught an error!')
            print(e)
        },
        warning = function(w){
            message('Caught an warning!')
            print(w)
        },
        finally = {
            message('All done, quitting.')
        }
    )
}

# calculate pi1 and return table
get_pi_1 <- function(df, sourceName, targetName){
    
    res <- run_qvalue(df$pvalue_target)

    # if P-values are all highly significant then qvalue will throw an error
    if("error" %in% class(res)){
        pi1 <- NA
    }
    if("pi0" %in% names(res) ){
        pi1 <- 1 - res$pi0
    }

    # create results table

    res_df <- data.frame( source_name = sourceName, target_name = targetName, pi1 = pi1)

    print(res_df)
    return(res_df)
}


library(optparse)
option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output folder', default = "."),
        make_option(c('-s', '--sourceName'), help = "the name of the source P-value distribution", default = "source"),
        make_option(c('-t', '--targetName'), help = 'the name of the target P-value distribution', default = "target"),
        make_option(c('-q', '--threshold'), help = 'how to threshold significant associations in source dataset', default = 0.05)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFolder <- opt$outFolder
sourceName <- opt$sourceName
targetName <- opt$targetName
threshold <- opt$threshold  

main <- function(){

#sourceName <- "Microglia_THA"
#targetName <- "Microglia_SVZ"
outFile <- file.path( outFolder, paste0(sourceName, ":", targetName, ":", threshold, ".qvalue.tsv") )

qtl_source <- pullData(sourceName, type = "QTL")
qtl_target <- pullData(targetName, type = "QTL")

# get source permutation QTLs (top QTLs)
source_top_res <- extractTopQTL(qtl_source,qvalue_threshold = threshold)

if( qtl_target$full_chr_type == "chr1" ){
    chroms <- paste0("chr", seq(1,22) ) 
}else{ 
    chroms <- 1:22
}

# iterate through chromosomes to get nominal results in target QTLs
target_res <- purrr::map_df( chroms, ~{extractTargetQTL(qtl_target, .x, source_top_res) })

# join together source and target
qvalue_res <- left_join(source_top_res, target_res, by = c("pheno", "snp") )

# remove missing values
qvalue_res <- qvalue_res[ !is.na(qvalue_res$pvalue_target), ]

# write data
save(qvalue_res, source_top_res, target_res,  file = gsub("tsv", "RData", outFile) )



# get pi1 of target compared to source
final_res <- get_pi_1( qvalue_res, sourceName, targetName)

# write table
message(" * writing to: ", outFile)
readr::write_tsv(final_res, path = outFile)
}

main()


# pull entry from GWAS/QTL DB
pullData <- function(dataset, type = "QTL"){
    message(Sys.time()," * selected dataset: ", dataset)
    db_path <- "/sc/arion/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"

    message(Sys.time()," * reading GWAS database from ", db_path)
    stopifnot( file.exists(db_path) )
    if( type == "GWAS"){
        n_sheet <- 3
    }
    if( type == "QTL"){
        n_sheet <- 2
    }
    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = n_sheet,na= c("", "-","NA")))

    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )

    #stopifnot( all( c("full_path", "full_chrom", "full_pos") %in% names(gwas) ) )

    #stopifnot( file.exists(gwas$full_path) )

    return(gwas)
}


# per chrom (or chunk?) read in QTL and extract SNP, phenotype, beta and SE
extractTargetQTL <- function(qtl, chr ){
    require(dplyr)
    stopifnot( file.exists(qtl$full_path) )

    # check whether tabixed with "1" or "chr1"
    if( qtl$full_chrom_type == "1.0" | qtl$full_chrom_type == 1  ){
        chr <- gsub("chr", "", chr)
    }else if( qtl$full_chrom_type == "chr1"){
        if( !grepl("chr", chr) ){
            chr <- paste0("chr", chr)
        }
    }

    # read in QTL
    # use dark awk magic to remove the clusterID while reading in
    if( qtl$phenotype == "sQTL" ){
        cmd <- paste0("ml bcftools; tabix ", qtl$full_path, " ", chr, ' | awk \'{gsub(/clu_[0-9]+_[+-]:/, ""); print }\' ')
    }else{
        cmd <- paste0("ml bcftools; tabix ", qtl$full_path, " ", chr)
    }
    message( " * ", cmd)
    result <- data.table::fread( cmd = cmd, nThread = 8 , fill = TRUE)
    message( " * QTL read in!" )
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
    names(result)[names(col_dict) == betaCol]   <- "beta"
    names(result)[names(col_dict) == seCol]   <- "se"
    
    # weird exception - Young Microglia stores p-values as log10 - convert
    if( qtl$full_p == "log10_p"){
        result$pvalue <- 10^result$pvalue
    }
    #return(result)

    # adjust EnsemblID - remove tag number
    if( qtl$phenotype == "eQTL" ){
        result$pheno <- stringr::str_split_fixed(result$pheno, "\\.", 2)[,1]
    }
    # merge SNP and pheno
    result$snp_gene <- paste0(result$snp ,"-", result$pheno)
    
    # subset out just SNP-Gene and P-value
    res_subset <- 
        result %>%
        select( snp_gene, pvalue)
    print( head(res_subset) )
    return(res_subset)
}


library(optparse)

option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = ""),
        make_option(c('-c', '--chromosome'), help = "which chromosome to read in", default = "")
)

option.parser <- OptionParser(usage = "%prog [options] QTL_dataset_name_1 QTL_dataset_name_2 ...QTL_dataset_name_N", description = "a script that combines QTL dataset into the format required by pisquared", option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = TRUE)
data_list <- opt$args
options <- opt$options

outFolder <- options$outFolder
chromo <- options$chromosome

if(chromo == ""){chromo <- NA}

message( " * chr is ", chromo )
message(" * outFolder is", outFolder) 

stopifnot( length( data_list) > 1 ) 

library(tidyverse)

# for each chromosome and each dataset, extract QTL beta and se, and merge by SNP-Gene pair
#data1 <- "Cervical_Spinal_Cord_NYGC"
#data2 <- "GTEX_Brain_Spinal_cord_cervical_c-1"
#data_list <- list(data1,data2) 


qtl_list <- purrr::map(data_list, pullData, type = "QTL" )

prefix <- paste( data_list, collapse = "-" )

outFile_final <- file.path(outFolder, paste0( prefix, ".metasoft.input.tsv" ) )
# if chromosome is specified then just read in from that chromosome
# else iterate through all

if( is.na(chromo) ){
    chr_nums <- 1:22
    chrs <- paste0("chr", chr_nums)
}else{
    chrs <- chromo
}

if( !dir.exists(outFolder) ){ dir.create(outFolder) }

for( chr in chrs){
    message(" * merging data in ", chr )
    outFile <- file.path( outFolder, paste0(prefix, ".", chr, ".pisquared.input.tsv.gz" ) )
    res <- purrr::map( qtl_list, ~{
        extractTargetQTL(.x, chr = chr)
    })  %>% 
    reduce( inner_join, by = "snp_gene"  )
    
    res <- res[ complete.cases(res),]     

    # give each column a unique name
    data_names <- unlist(purrr::map( data_list, ~{ res <- paste( c("pvalue"), .x, sep = "."); return(res) }))
    colnames(res) <- c("snp_gene", data_names)
    
    print(head(res) )
    print( names(res) )
    message( " * writing to ", outFile )  
    if( as.character(chr) == "1"){
        write_tsv( res, path = outFile, col_names = TRUE)
    }else{
        write_tsv( res, path = outFile, col_names = FALSE)
    }
    rm(res)
    gc()
}



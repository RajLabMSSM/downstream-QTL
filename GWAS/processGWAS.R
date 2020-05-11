
# Process GWAS
# Jack Humphrey

## For a GWAS summary stat table
## Split by chromosome
## Sort by position
## Write out as TSV
## bgzip and tabix index

## WARNING: this script does not lift over coordinates or flip alleles. You should have done this already.

library(data.table)
setDTthreads(0) # set number of threads to number provided
library(readr)
library(dplyr)
library(optparse)

# REFACTORED TO EXCLUSIVELY USE THE RAJ LAB GWAS DATABASE

pullGWAS <- function(dataset){
    message(Sys.time()," * selected dataset: ", dataset)
    db_path <- "/sc/hydra/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"
    
    message(Sys.time()," * reading GWAS database from ", db_path)
    stopifnot( file.exists(db_path) )

    gwas_db <- suppressMessages(readxl::read_excel(db_path, sheet = 2))

    stopifnot( dataset %in% gwas_db$dataset )
    
    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )
    
    stopifnot( all( c("full_path", "full_chrom", "full_pos") %in% names(gwas) ) )

    stopifnot( file.exists(gwas$full_path) )
    
    return(gwas) 
}

# read in header of summary stat file
# get out mapping of column names to column numbers
columnDictionary <- function(file_path){
  stopifnot(file.exists(file_path) )
  # Get the index of each column name
  f <- data.table::fread(file_path, nrows = 0, header = TRUE)
  cNames <- colnames(f)
  colDict <- setNames(1:length(cNames), cNames  )
  return(colDict)
}

# read in summary stats
# if "chr" present in col_chr then remove it
# feature for later: if no chr or pos column present then match on RS ID
# sort by chr and position
# write out sorted file to reference folder
# tabix index the file
sortGWAS <- function(gwas, out_folder = "./"){
    dataset <- gwas$dataset
    
    col_chr <- gwas$full_chrom
    col_pos <- gwas$full_pos    
    message(Sys.time()," * reading in GWAS file")
    gwas_df <- data.table::fread(gwas$full_path, nThread = 4)

    # sample chr column - if "chr" present then remove
    gwas_df[[col_chr]] <- gsub("chr", "", gwas_df[[col_chr]] )

    gwas_df[[col_chr]] <- as.numeric(gwas_df[[col_chr]] )
    gwas_df[[col_pos]] <- as.numeric(gwas_df[[col_pos]] )
    # sort GWAS by chr and pos
    message(Sys.time()," * sorting GWAS")
    gwas_sorted <- gwas_df[order(get(col_chr), get(col_pos) )] # get() makes R evaluate the column name as a variable not a string, to allow sorting

    out_path <- paste0(out_folder, dataset, ".processed.tsv")
    message(Sys.time()," * writing sorted GWAS file to ", out_path)
    readr::write_tsv(gwas_sorted, out_path)
    stopifnot(file.exists(out_path) )
}

tabixGWAS <- function(gwas, outFolder = "./"){
    dataset <- gwas$dataset

    out_path <- paste0(outFolder, dataset, ".processed.tsv")
    stopifnot(file.exists(out_path) )
     
    col_chr <- gwas$full_chrom
    col_pos <- gwas$full_pos

    col_dict <- columnDictionary(gwas$full_path)

    # get column numbers
    n_chr <- which(names(col_dict) == col_chr )
    n_pos <- which(names(col_dict) == col_pos )
    
    # bgzip
    message(Sys.time()," * bgzipping GWAS" )
    system( paste("ml bcftools; bgzip -f ", out_path) )

    out_path_gzip <- paste0(out_path, ".gz")
    stopifnot(file.exists(out_path_gzip) )

    # tabix file
    message(Sys.time()," * tabixing GWAS" )
    tabix_cmd <- paste0("ml bcftools; tabix -f -S 1 -s ", n_chr, " -b ", n_pos, " -e ", n_pos, " ", out_path_gzip)
    system(tabix_cmd)
    out_path_tabix <- paste0(out_path_gzip, ".tbi")
    stopifnot( file.exists( out_path_tabix) )
    message(Sys.time()," * processing complete!" )
}


option_list <- list(
    make_option(c('-n', '--dataset' ), help='the name of the dataset. This must match the value in the database', default = "Marioni_test"),
    make_option(c('-o', '--out_folder'), help = "the full path to where the processed GWAS should be written to", default = "/sc/hydra/projects/ad-omics/data/references/GWAS/")
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

dataset <- opt$dataset
out_folder <- opt$out_folder


message(Sys.time()," * GWAS summary stat processor ")

gwas <- pullGWAS(dataset)

sortGWAS(gwas, out_folder)
tabixGWAS(gwas, out_folder)





# Process GWAS
# Jack Humphrey

## For a GWAS summary stat table
## Split by chromosome
## Sort by position
## Write out as TSV
## bgzip and tabix index

## WARNING: this script does not lift over coordinates or flip alleles. You should have done this already.

library(data.table)
library(readr)
library(dplyr)
library(optparse)

#option_list <- list(
#    make_option(c('-i', '--inFile'), help = "the full GWAS summary stats"),
#    make_option(c('-o', '--outFile'), help='the path to the output file - without file type suffix', default = ""),
#    make_option(c('--chrCol'), help = "the column number that stores the chromosome of the variant", default = 3),
#    make_option(c('--posCol'), help = "the column number that stores the genomic position of the variant", default = 4),
#    make_option(c('--noChrPrefix'), help = "don't prepend chr to the chromosome column values", action="store_true", default=FALSE) 
#)

# REFACTORED TO EXCLUSIVELY USE THE RAJ LAB GWAS DATABASE

option_list <- list(
    make_option(c('--dataset' ), help='the name of the dataset. This must match the value in the database', default = "")
)


option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

dataset <- opt$dataset



pullGWAS <- function(dataset){

    db_path <- "/sc/hydra/projects/ad-omics/data/references/GWAS/GWAS-QTL_data_dictionary.xlsx"
    stopifnot( file.exists(db_path) )

    gwas_db <- readxl::read_excel(db_path, sheet = 2)

    stopifnot( dataset %in% gwas_db$dataset )

    gwas <- gwas_db[ gwas_db$dataset == dataset & !is.na(gwas_db$dataset), ]
    stopifnot( nrow(gwas) == 1 )
    
    stopifnot( all( c("full_path", "full_chrom", "full_pos") %in% names(gwas) ) )
    stopifnot( all( !is.na( c(gwas$full_path, gwas$full_chrom, gwas$full_pos ) ) ) )
    
    full_chrom <- gwas$full_chrom
    full_pos <- gwas$full_pos
    full_path <- gwas$full_path

    stopifnot( file.exists(full_path) )
    
    out <- list(full_chrom = full_chrom, full_pos = full_pos, full_path = full_path)
    return(out) 
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

print(" * GWS summary stat processor ")

if( noChrPrefix == TRUE){
print(" * no adding chr string to chromosomes")
}else{
print(" * prepending \"chr\" to chromosome names")
}

print(" * reading in file")

gwas <- pullGWAS(dataset)

col_dict <- columnDictionary(gwas$full_path)

#gwas_df <- data.table(fread(gwass$full_path, threads = 4))

col_chr <- gwas$full_chrom
col_pos <- gwas$full_pos

# get column numbers
n_chr <- which(names(col_dict) == col_chr )
n_pos <- which(names(col_dict) == col_pos )

# chr prefix - insist on it?

if( noChrPrefix == FALSE ){
gwas_df[[col_chr]] <- paste0("chr", gwas_df[[col_chr]] )
}

names(gwas_df)[num_chr] <- "chr"
names(gwas_df)[num_pos] <- "pos"

setkey(gwas_df, "chr")
# get unique values of chromosome column
all_chrs <- unique( gwas_df[,chr] ) 

# CHR NAME ISSUE

sortFile <- function(gwas, col_dict){
    n_chr <- which(names(col_dict) == gwas$full_chrom )
    n_pos <- which(names(col_dict) == gwas$full_pos )
}


tabixFile <- function(gwas, col_dict){
    n_chr <- which(names(col_dict) == gwas$full_chrom )
    n_pos <- which(names(col_dict) == gwas$full_pos )
    
    
}


for(i in all_chrs){
    print(paste0(" * processing ", i ))
    chr_df <- gwas_df[i]
    # sort by position
    chr_df_sorted <- chr_df[order(pos)]
    
    # remove duplicate entries - seem to exist in EBI GWAS
    chr_df_final <- chr_df_sorted %>% as.data.frame() %>% distinct()
    
    # write out file
    out_path <- paste0(outFile, "_", i, ".tsv")

    write_tsv( chr_df_final, path = out_path)
    
    # bgzip file
    system( paste("ml bcftools; bgzip ", out_path) )

    out_path_gzip <- paste0(out_path, ".gz")
    stopifnot(file.exists(out_path_gzip) )

    # tabix file
    tabix_cmd <- paste0("ml bcftools; tabix -f -S 1 -s ", num_chr, " -b ", num_pos, " -e ", num_pos, " ", out_path_gzip) 
    system(tabix_cmd)
    out_path_tabix <- paste0(out_path_gzip, ".tbi")
    stopifnot( file.exists( out_path_tabix) )
}

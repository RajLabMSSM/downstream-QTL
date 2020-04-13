
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

option_list <- list(
    make_option(c('-i', '--inFile'), help = "the full GWAS summary stats"),
    make_option(c('-o', '--outFile'), help='the path to the output file - without file type suffix', default = ""),
    make_option(c('--chrCol'), help = "the column number that stores the chromosome of the variant", default = 3),
    make_option(c('--posCol'), help = "the column number that stores the genomic position of the variant", default = 4),
    make_option(c('--noChrPrefix'), help = "don't prepend chr to the chromosome column values", action="store_true", default=FALSE) 
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
inFile <- opt$inFile
chrCol <- opt$chrCol
posCol <- opt$posCol
noChrPrefix <- opt$noChrPrefix

# read in full GWAS 
#gwas <- "../example/nicolas_als_chr21_gwas.tsv.gz"
#out_prefix <- "nicolas_als"

if(!is.na(chrCol) & !is.na(posCol) ){
    num_chr <- chrCol
    num_pos <- posCol
}

print(" * GWS summary stat processor ")

if( noChrPrefix == TRUE){
print(" * no adding chr string to chromosomes")
}else{
print(" * prepending \"chr\" to chromosome names")
}

print(" * reading in file")
gwas_df <- fread(inFile)

col_chr <- names(gwas_df)[num_chr]
col_pos <- names(gwas_df)[num_pos]

if( noChrPrefix == FALSE ){
gwas_df[[col_chr]] <- paste0("chr", gwas_df[[col_chr]] )
}

names(gwas_df)[num_chr] <- "chr"
names(gwas_df)[num_pos] <- "pos"

setkey(gwas_df, "chr")
# get unique values of chromosome column
all_chrs <- unique( gwas_df[,chr] ) 

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

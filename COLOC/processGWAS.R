library(data.table)
library(readr)
library(dplyr)
library(optparse)

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file - without file type suffix', default = ""),
    make_option(c('-i', '--inFile'), help = "the full GWAS summary stats"),
    make_option(c('-t', '--type'), help = "where the GWAS is from. These can be \"EBI\" - more to be added", default = "EBI")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
inFile <- opt$inFile
type <- opt$type


# read in full GWAS 
#gwas <- "../example/nicolas_als_chr21_gwas.tsv.gz"
#out_prefix <- "nicolas_als"

if( type == "EBI"){
    col_chr <- "hm_chrom"
    col_pos <- "hm_pos"
    num_chr <- 3
    num_pos <- 4
    chr_append <- TRUE # should "chr" be appended to the chromosome values?
}


print(" * GWS summary stat processor ")
print(" * reading in file")
gwas_df <- fread(inFile)

if( chr_append == TRUE ){
gwas_df[[col_chr]] <- paste0("chr", gwas_df[[col_chr]] )
}

names(gwas_df)[which(names(gwas_df) == col_chr )] <- "chr"
names(gwas_df)[which(names(gwas_df) == col_pos )] <- "pos"

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

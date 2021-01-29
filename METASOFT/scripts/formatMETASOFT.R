library(data.table)
library(optparse)

option_list <- list(
    make_option(c('--vcf', '-v' ), help='The full path to the VCF used in the QTL analysis', default = ""),
    make_option(c('-o', '--output'), help = "the full path to the output file", default = ""),
    make_option(c('-i', '--input'), help = "the input metasoft file", default = ""),
    make_option(c('-c', '--chromo'), help = "the chromosome used", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

vcf <- opt$vcf
input <- opt$input
output <- opt$output
chromo <- as.numeric(opt$chromo)

# for testing
#input <- "/sc/hydra/projects/ad-omics/microglia_omics/METASOFT/sQTL/Microglia_sQTL_THA-Microglia_sQTL_SVZ-Microglia_sQTL_MFG-Microglia_sQTL_STG.13.metasoft.res.tsv"
#output <- "test.tsv"
#vcf <- "/sc/hydra/projects/pd-omics/glia_omics/eQTL/post_imputation_filtering/eur/filtered_variants/AllChr.hg38.sort.filt.dbsnp.snpeff.vcf.gz"
#vcf <- "/sc/hydra/projects/ad-omics/microglia_omics/METASOFT/sQTL/vcf.processed.txt"
#chromo <- 13

stopifnot(!is.na(chromo) )


#### format METASOFT output

# read in metasoft file
message(" * reading in METASOFT results")
df <- data.table::fread( input, nThread = 8)
# split SNP and Phenotype

# RE2C output calls SNPs rsID, RECOV calls them "snp"
names(df)[1] <- "RSID"
df[, c("variant_id", "phenotype_id") := tstrsplit(RSID, "-", fixed=TRUE)]
# set key as variant_id
setkey(df, variant_id)

# use VCF to match SNP to coordinate
stopifnot(file.exists(vcf) )

# read in SNP positions from VCF

message(" * reading in SNP positions from VCF")
snp_df <- data.table::fread(vcf, col.names = c("chr", "pos", "variant_id"), key = "variant_id", nThread = 8)

message(" * merging files on shared RSID" )
# inner join between the VCF SNPs and the meta-analysis
all_df <- merge(snp_df,df)

rm(snp_df)
rm(df)

# sort by position
setkey(all_df, pos)
#all_df <- all_df[sort(pos)]

if( chromo == 1){ 
    write_cols <- TRUE
}else{ 
    write_cols <- FALSE 
}

message(" * saving file to ", output)

# write out using data.table fwrite
fwrite( all_df, file = output, sep = "\t", na = NA, quote = FALSE, col.names = write_cols, nThread = 8 )


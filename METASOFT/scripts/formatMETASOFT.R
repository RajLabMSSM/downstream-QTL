library(data.table)
library(optparse)

option_list <- list(
    make_option(c('--vcf', '-v' ), help='The full path to the VCF used in the QTL analysis', default = ""),
    make_option(c('-o', '--output'), help = "the full path to the output prefix", default = ""),
    make_option(c('-i', '--input'), help = "the input metasoft file", default = ""),
    make_option(c('-c', '--chromo'), help = "the chromosome used", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

vcf <- opt$vcf
input <- opt$input
prefix <- opt$output
chromo <- as.numeric(opt$chromo)

# for testing
#input <- "/sc/hydra/projects/ad-omics/microglia_omics/METASOFT/sQTL/Microglia_sQTL_THA-Microglia_sQTL_SVZ-Microglia_sQTL_MFG-Microglia_sQTL_STG.13.metasoft.res.tsv"
#output <- "test.tsv"
#vcf <- "/sc/hydra/projects/pd-omics/glia_omics/eQTL/post_imputation_filtering/eur/filtered_variants/AllChr.hg38.sort.filt.dbsnp.snpeff.vcf.gz"
#vcf <- "/sc/hydra/projects/ad-omics/microglia_omics/METASOFT/sQTL/vcf.processed.txt"
#chromo <- 13

stopifnot(!is.na(chromo) )

output_full <- paste0(prefix, ".full_assoc.tsv")
output_top <- paste0(prefix, ".top_assoc.tsv")

#### format METASOFT output

# read in metasoft file
message(" * reading in METASOFT results")
df <- data.table::fread( input, nThread = 8)
#cripts/formatMETASOFT.Rsplit SNP and Phenotype

# RE2C output calls SNPs rsID, RECOV calls them "snp"
names(df)[1] <- "RSID"
df[, c("variant_id", "phenotype_id") := tstrsplit(RSID, "-", fixed=TRUE)]

## multiple testing correction
# For METASOFT, correct fixed and random effect Ps per gene by BH
# for each phenotype calculate and add adjusted P values
if( "PVALUE_RE2" %in% names(df) ){
    df[, c("PADJ_RE2", "PADJ_FE") := .(p.adjust(.SD$PVALUE_RE2, method = "bonferroni"), p.adjust(.SD$PVALUE_FE, method = 'bonferroni') ), by = phenotype_id]
}
if( "RE2Cp" %in% names(df) ){
    df[, "PADJ_RE2C" := p.adjust(.SD$RE2Cp, method = "bonferroni") , by = phenotype_id]

}
# output top association per gene

## data.table - iterate through groups in data frame using x[, print(.SD), by = phenotype_id]


# get top assoc x[order(PVALUE_RE), head(.SD,1), by = phenotype_id]


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

message(" * saving file to ", prefix)

# get top assocations
if( "PADJ_RE2" %in% names(all_df) ){
    top <- all_df[order(PADJ_RE2), head(.SD,1), by = phenotype_id]
    }
# for RE2C
if( "RE2Cp" %in% names(all_df) ){
     top <- all_df[order(PADJ_RE2C), head(.SD,1), by = phenotype_id]
}
#save.image("debug.RData")
# write out using data.table fwrite
fwrite( top, file = output_top, sep = "\t", na = NA, quote = FALSE, col.names = write_cols, nThread = 1)

fwrite( all_df, file = output_full, sep = "\t", na = NA, quote = FALSE, col.names = write_cols, nThread = 8 )


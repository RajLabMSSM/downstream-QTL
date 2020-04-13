
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

gwas <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/processed/Nicolas_2018_processed_"
hits_file <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/COLOC/Nicolas_2018/Nicolas_2018_hits_1e-7.tsv"

hits <- read_tsv(hits_file)

makeCoords <- function(df, flank = 1e6){
    stopifnot( all(c("chr", "pos") %in% names(df) ) )
    stopifnot( flank > 0)
    coords <- paste0(df$chr, ":", df$pos - flank, "-",df$pos + flank) 
    return(coords)
}

extractGWAS <- function(coords, gwas, chr){
    # assume coord is a string following chr:start-end format
    cmd <- paste("ml bcftools; tabix ", gwas, chr, coords ) 
    result <- system(cmd, intern = TRUE)
    return(result)
}





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


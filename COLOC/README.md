# COLOC scripts

## Inputs:

* QTL nominal summary statistics (per-chromosome parquet files)

Created by TensorQTL

* GWAS summary statistics

Use EBI GWAS Catalogue as specification


### Functions:

**indexGWAS** - tabix index a GWAS summary stat file

    input: GWAS file

    output: tabix indexed GWAS file

    parameters: which columns to use - use EBI GWAS catalogue as default

**extractTopHits** - from a summary stat file, get out the top hit variants at a given significance threshold (default 1e-8)

    input: tabix indeed GWAS file

    output: table with coordinates, betas and P-values of top hits

    parameters: given significance threshold


**extractGWAS** - from a summary stat file, extract GWAS statistics within a coordinate range

    input: tabix indexed GWAS file

    output: list of P-values and coordinates

    parameters: which column is P-value (default to EBI GWAS catalogue)


**extractNominalQTLs** - from a set of parquet files, extract QTL statistics within a coordinate range

    input: folder containing per-chromosome parquet files containing nominal QTL results

    output: list of P-values and coordinates

    parameters: which column is p-value (default should be the same except when using interaction results)


**runCOLOC** - run COLOC on with the two P-value distributions

    input: 2 lists of P-values and coords, one from GWAS and one from QTL

    output: summary of COLOC results

    parameters: ???


## Putting steps together

0. Index GWAS

1. Extract top hits from GWAS

2. For a given GWAS and a top hit variant:
    - extract all nominal QTL associations within 1MB from the hit variant (ExtractNominalQTLs)
    - extract all GWAS stats within the same region
    for each gene:
        run COLOC to test for colocalisation with each gene's QTL P-value distribution and the GWAS distribution



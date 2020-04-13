# COLOC scripts

## Inputs:

* QTL nominal summary statistics (per-chromosome parquet files)

Created by TensorQTL

* GWAS summary statistics

Use EBI GWAS Catalogue as specification

## Dependencies:

### Utilities:

tabix (ml bcftools/1.9)

### R packages

* arrow 

```
git clone https://github.com/apache/arrow
cd arrow/r
R CMD INSTALL 
```

* COLOC

### Functions:

**processGWAS.R** - take a GWAS summary stat file and split by chromosome, sort by position and tabix index

    input: GWAS file

    output: tabix indexed GWAS files

    parameters: which columns to use - use EBI GWAS catalogue as default

```
Options:
    -i INFILE, --inFile=INFILE
        the full GWAS summary stats

    -o OUTFILE, --outFile=OUTFILE
        the path to the output file - without file type suffix

    --chrCol=CHRCOL
        the column number that stores the chromosome of the variant

    --posCol=POSCOL
        the column number that stores the genomic position of the variant

    --noChrPrefix
        don't prepend chr to the chromosome column values

    -h, --help
        Show this help message and exit
```

This will create a set of individual tabixed chromosome files.

**extractTopHits.R** - from a summary stat file, get out the top hit variants at a given significance threshold (default 1e-8)

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



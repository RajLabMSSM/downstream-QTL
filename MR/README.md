# MR scripts

## Inputs:

* QTL nominal summary statistics (per-chromosome parquet files)

Created by TensorQTL

* GWAS summary statistics

Use EBI GWAS Catalogue as specification

## Dependencies:

### Utilities:

tabix (ml bcftools/1.9)

### R packages


* MR

### Functions:


**run_MR.R** - run MR on with the two P-value and beta distributions

```
Usage: Rscript scripts/run_MR.R [options]


Options:
    -o OUTFOLDER, --outFolder=OUTFOLDER
        the path to the output file

    -g GWAS, --gwas=GWAS
        the dataset ID for a GWAS in the GWAS/QTL database

    -q QTL, --qtl=QTL
        the dataset ID for a QTL dataset in the GWAS/QTL database

    -f FDR, --fdr THREHSOLD
	the exposure filter option, defalt=FALSE, which use GWAS threshold 5e-8, if TRUE, use FDR < 0.05

    --debug
        load all files and the nsave RData without running MR

    -h, --help
        Show this help message and exit
```


## Putting steps together

0. Index GWAS

1. Extract top hits from GWAS

These are now taken care of in the GWAS scripts

2. For a given GWAS and a top hit variant:
    - extract all nominal QTL associations within 1MB from the hit variant (ExtractNominalQTLs)
    - extract all GWAS stats within the same region
    - filter by threshold GWAS threshold or FDR
    for each gene:
        run MR to test for colocalisation with each gene's QTL P-value distribution and Beta and the GWAS distribution


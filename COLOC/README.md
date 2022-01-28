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


* COLOC

### Functions:


**runCOLOC.R** - run COLOC on with the two P-value distributions

```
Usage: Rscript scripts/run_COLOC.R [options]


Options:
    -o OUTFOLDER, --outFolder=OUTFOLDER
        the path to the output file

    -g GWAS, --gwas=GWAS
        the dataset ID for a GWAS in the GWAS/QTL database

    -q QTL, --qtl=QTL
        the dataset ID for a QTL dataset in the GWAS/QTL database

    --debug
        load all files and the nsave RData without running COLOC

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
    for each gene:
        run COLOC to test for colocalisation with each gene's QTL P-value distribution and the GWAS distribution


## LDlinkR

to calculate LD between the lead GWAS and lead QTL SNPs, we use the LDLinkR R package. You need to generate an API access token and add it to your ~/.Renviron
`LDLINK_TOKEN=" "`
[see here](https://cran.r-project.org/web/packages/LDlinkR/vignettes/LDlinkR.html)


# Working with GWAS

## Rationale

We are attempting to maintain a database of all useful GWAS for the lab. This is mostly Neurological and Neuropsych disorders, plus some immune and other traits.

## Methodology

The idea is that we download the raw summary stats to Minerva. We then enter the details of the GWAS and the file to the database.

Using the database we can then preprocess the GWAS summary statistics. Right now this involves sorting and tabix indexing for random access, but could expand to lift over, allele flipping etc.

First enter your GWAS summary stats as an entry to the [GWAS-QTL database](https://drive.google.com/file/d/1BgLQaRZd9L7JoO8IbpzhUCRFdTdLgJvD/view?usp=sharing). 

The full path to the raw summary stats should be entered in the **full_path** column.

Give the GWAS a unique name (eg Lambert_2013).

Fill in the spreadsheet with the metadata about the GWAS and the names of the columns:

### GWAS metadata

dataset: Name of the dataset
phenotype:   Phenotype being associated (eg Alzheimer's Disease)
build:   Human genome build (eg hg19 or hg38)
type:    Whether case/control "cc" or quantitative trait "quant" design used
N:   Number of individual donors in analysis
prop_cases:  Proportion of donors that are cases (in case/control only)
top_path:    Full path to top associations. If you don't have this we can create it.
top_filtered_path:   Full path to top associations  with loci filtered out - see notes
top_sheet:   If top associations are in an Excel workbook, then which worksheet should we use? Either name of sheet of the sheet number
top_path:    Absolute path to top associations / summary stats
full_path:   Absolute path to full associations / summary stats
full_processed_path: Absolute path to full associations / summary stats; after being tabix-indexed, having missing RSIDs filled, and inferring SE.
reference:   URL to paper, if applicable

### Column data:

full_snp:     Name of the SNP RSID column in full summary associations
full_chrom:  Name of chromosome column in full summary associations
full_pos:    Name of position column in full summary associations
full_effect: Name of the effect size column in full summary associations.  While beta is often the reported effect size, other metrics may be reported in addition or instead. 
full_se: Name of standard error of the effect size column in full summary associations
full_p:  Name of P-value column in full summary associations
full_stat:   Name of the t-statistic column in the full summary associations. This can be used to infer standard error (effect / t-stat = SE)."
full_freq:   The name of the column with the frequency of the effect allele. Can be used to infer MAF.
full_maf:    Name of Minor Allele Frequency (Effect allele frequency) in full summary association
full_A1: The allele that the valence of the effect size is relative to.
full_A2: The other allele.
full_alleleTypes:   What full_A1 and full_A2 correspond to.


Once you've entered your GWAS, make sure that you synchronise the GWAS/QTL database to the cluster - ask Jack!

Jack - at some point work out how to access remote Google Sheet from minerva.

For now, sync_db.sh will copy over your local version of the database to its place on minerva and grant group read permission


## Harmonising summary stats

```
ml R/3.6.0
Rscript processGWAS.R -n <dataset> -o <outFolder> 

```

<dataset> must match the dataset name you put in the GWAS/QTL database.

Currently I'm keeping the processed GWAS sumstats here: /sc/arion/projects/ad-omics/data/references/GWAS

## Getting top loci

TO BE WRITTEN


## Applications of processed GWAS

Again, using the database we can then perform various tasks with the GWAS summary stats, including:

* Fine-Mapping with echolocatoR

* Colocalisation analysis with COLOC

* LD-score regression, mediation analysis, TWAS and so much more.


## Summary

Indexing a GWAS in the database and pre-processing it is the start of a beautiful and well-organised science adventure for the whole lab to enjoy!

[Link to GWAS database](https://drive.google.com/file/d/1BgLQaRZd9L7JoO8IbpzhUCRFdTdLgJvD/view?usp=sharing)

To do: how to process GWAS summary stats

# qvalue QTL sharing

Applying the π0 method to assess rates of sharing between two QTL results

Inputs:

* Source QTL dataset
    - Permutation results
* Target QTL dataset
    - Full nominal results

Algorithm:

For source dataset, extract SNP-Gene pairs for all significant genes according to the permutation result - parse a text file

For target dataset, extract those same SNP-Gene pairs from the nominal results (parse the parquet files)

Join the two tables together by shared SNP-Gene pair columns

Run qvalue on the target P-value column to get Pi1. Take 1 - Pi1 to get Pi0.

## Running

This interfaces with the [Raj Lab GWAS/QTL database](https://docs.google.com/spreadsheets/d/1BgLQaRZd9L7JoO8IbpzhUCRFdTdLgJvD/edit#gid=473131646) - this is where we keep a standard set of metadata for all GWAS and QTL summary stats that we have access to.

The qvalue pi1 can be calculated between any two QTL summary statistics available, with the caveat that only QTL datasets that have permutation results can be used as the "source" set.

For full list of datasets, check [here](https://docs.google.com/spreadsheets/d/1BgLQaRZd9L7JoO8IbpzhUCRFdTdLgJvD/edit#gid=473131646)

```
Rscript run_qvalue --source <name of source dataset> --target <name of target dataset> --outFolder <path to store results>
```

Merging multiple results can be done with 

```
merge_qvalue.R --outFolder <path where results are>
```

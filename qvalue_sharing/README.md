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

currently just uses tensorQTL outputs - I will expand to use other sources of summary stat

This assumes that permutation results are kept in a .cis_qtl.txt.gz file, and the nominal p values are kept in parquet files, one per chromosome.

--source refers to the permutation p-values

--target refers to the nominal p-values

```
# example using NYGC ALS QTLs
qtl=/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/

# first prepare data
python prepare_qvalue_sharing.py --source $qtl/FrontalCortex_expression/peer30/FrontalCortex_expression_peer30.cis_qtl.txt.gz --target $qtl/Cerebellum_expression/peer30/ -o test.tsv

# then run qvalue on output file
Rscript run_qvalue.R --inFile test.tsv --outFile qvalue_res.txt --sourceName FrontalCortex --targetName Cerebellum

```

qvalue.res.txt is a single row dataFrame containing the values "source", "target" and "pi1"

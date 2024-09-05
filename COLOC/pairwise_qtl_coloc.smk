import sys
import pandas as pd
import glob
import re
import os

# Pairwise QTL COLOC pipeline
# Jack Humphrey
outFolder = config["outFolder"]
inFolder = config["inFolder"]

target_files = [ x for x in glob.glob(inFolder + "*genewide_QTL_pairs*tsv.gz" )]
target_ids = [os.path.basename(re.sub('.tsv.gz', '', x, count=1)) for x in target_files ]


shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge; conda activate snakemake;")
rule all:
    input:
        expand( outFolder + "{target}_COLOC.tsv.gz", target = target_ids)

rule pairwise_QTL_COLOC:
    input:
        inFolder + "{target}.tsv.gz"
    output:
        outFolder + "{target}_COLOC.tsv.gz"
    params:
        script = "scripts/run_COLOC.R"
    shell:
        "ml R;" 
        "Rscript {params.script} --targets {input} -o {outFolder} --lowmem --threads 1"

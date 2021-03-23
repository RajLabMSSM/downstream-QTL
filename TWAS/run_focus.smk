# Run FOCUS to fine-map TWAS models across multiple panels

import os
import pandas as pd

data_code = config["data_code"]
sumstat_file = config["sumstat_file"]

fusion_dir = "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/"
out_folder = "results/focus/"
prefix = out_folder + data_code

metadata = pd.read_csv(config["weight_files"])
metadata_dict = metadata.set_index('dataset').T.to_dict()
weight_codes = metadata["dataset"]

#print(weight_codes)

ld_ref = fusion_dir + "LDREF_hg38"
ld_prefix = "hwe1e6.1000G.EURn404.GRCh38_fk.chr"

shell.prefix('export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge; conda activate focus; ml R/3.6.0;')



rule all:
    input:  out_folder + data_code + ".db"

rule FOCUS_create:
    input:
         metadata["pos_file"]
        # all pos files
    output:
        # focus database
        out_folder + data_code + ".db"
    params:
        outname = out_folder + data_code
    run:
        #shell("focus -h")
        for i in range(len(input)):
            tissue = metadata["dataset"][i]
            pos_file = metadata["pos_file"][i]
            
            shell("focus import --verbose  --tissue {tissue} --name {tissue} --assay rnaseq --output {params.outname} {pos_file} fusion")

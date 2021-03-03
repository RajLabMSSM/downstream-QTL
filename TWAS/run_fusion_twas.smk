import os
import pandas as pd

data_code = config["data_code"]
sumstat_file = config["sumstat_file"]

fusion_dir = "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/"
out_folder = "results/"
prefix = out_folder + data_code
sumstat_file = "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/nicolas_hg38_gwas_munged.sumstats.gz"

#weight_files = ["GTEx.Whole_Blood", "CMC.BRAIN.RNASEQ", "CMC.BRAIN.RNASEQ_SPLICING", "Brain_Cerebellum", "Brain_Spinal_cord_cervical_c-1"]

chromosomes = [i for i in range(1,23)]
# for testing
#chromosomes = [21]

metadata = pd.read_csv(config["weight_files"])
metadata_dict = metadata.set_index('dataset').T.to_dict() 
weight_codes = metadata["dataset"]
print(weight_codes)

ld_ref = fusion_dir + "LDREF_hg38"
ld_prefix = "hwe1e6.1000G.EURn404.GRCh38_fk.chr"



rule all:
    input: 
        expand( out_folder + data_code + ".{weight}.fusion_res.tsv", weight = weight_codes )

rule fusion_assoc_test:
    output:
         out_folder + "{weight}/" + data_code + ".{weight}.{chromo}.twas.txt"
    run:
        weight_dir = metadata_dict[wildcards.weight]["weight_dir"]
        pos_file = metadata_dict[wildcards.weight]["pos_file"]
        shell("ml R/3.6.0; \
        Rscript {fusion_dir}/FUSION.assoc_test.R --sumstats {sumstat_file}  \
        --weights {pos_file}  \
        --chr {wildcards.chromo}  \
        --weights_dir {weight_dir}  \
        --ref_ld_chr {ld_ref}/{ld_prefix}  \
        --out {output} ")
 
rule merge_fusion_res:
    input:
        expand( out_folder + "{weight}/" +  data_code + ".{weight}.{chromo}.twas.txt", chromo = chromosomes, allow_missing=True )
    output:
        out_folder + data_code + ".{weight}.fusion_res.tsv"
    
    params:
        script = "scripts/merge_fusion_res.R"
    
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {out_folder}{wildcards.weight}/ -o {output} "

# some datasets are missing pos files
# create them if missing by using GENCODE metadata and the Ensembl IDs used by the weights

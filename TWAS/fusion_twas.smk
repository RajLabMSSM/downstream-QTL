import os
import pandas as pd

data_code = config["data_code"]
weight_files = config["weight_files"]
sumstat_file = config["sumstats_file"]

fusion_dir = "/sc/arion/projects/ad-omics/data/software/fusion_twas-master"

out_folder = "results/"
prefix = out_folder + data_code
sumstat_file = "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/nicolas_hg38_gwas_munged.sumstats.gz"

weight_files = ["GTEx.Whole_Blood", "CMC.BRAIN.RNASEQ", "CMC.BRAIN.RNASEQ_SPLICING", "Brain_Cerebellum", "Brain_Spinal_cord_cervical_c-1"]

#chromosomes = [i for i in range(1,23)]
# for testing
chromosomes = [21]

rule all:
    input: 
        expand( out_folder + "/{weight}/" + data_code + ".{weight}.{chromo}.twas.txt", chromo = chromosomes, weight = weight_files )

rule fusion_assoc_test:
    output:
         out_folder + "/{weight}/" + data_code + ".{weight}.{chromo}.twas.txt"
    shell:
        "ml R/3.6.0;"
        "Rscript {fusion_dir}/FUSION.assoc_test.R --sumstats {sumstat_file} "
        "--weights {fusion_dir}/WEIGHTS/{wildcards.weight}/{wildcards.weight}.pos "
        "--chr {wildcards.chromo} "
        "--weights_dir {fusion_dir}/WEIGHTS/ "
        "--ref_ld_chr {fusion_dir}/LDREF/1000G.EUR. "
        "--out {output} "
 

import os
import pandas as pd


# for a QTL dataset, create FUSION weights for TWAS

# inputs (variables for each dataset)
# VCF - the full genotype matrix for common variants
# phenotype_df - this is produced by the QTL-mapping pipeline for expression or splicing phenotypes
# covariate_df - known covariates file - produced by QTL-mapping pipeline transposed to each covariate is a separate column
# data_code - the name of the dataset

# test data - the thoracic spinal cord eQTLs
#data_code = "NYGC_Thoracic_Spinal_Cord"
#vcf = "/sc/arion/projects/als-omics/WGS_QC/NYGC_Freeze02_European_Feb2020/WGS_QC_Pipeline/NYGC_Freeze02_European_Feb2020/output/chrAll_QCFinished_MAF0.01.anno.vcf.gz"
#mode = "eQTL"
#n_sample = 68
#pheno_df = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/ThoracicSpinalCord_expression/ThoracicSpinalCord_expression.phenotype.tensorQTL.gene.bed.gz"
#pheno_df = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/TWAS/nygc_cervical_test_pheno.bed"
#cov_df = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/ThoracicSpinalCord_expression/peer10/ThoracicSpinalCord_expression_peer10.gene.combined_covariates.txt"
#out_folder = "results/"


data_code = config["data_code"]
vcf = config["vcf"]
mode = config["mode"]
n_sample = config["n_sample"]
pheno_df = config["pheno_df"]
cov_df = config["cov_df"]
out_folder = config["out_folder"]


# hardcoded inputs
# LD reference - from FUSION
# number of chunks - how many phenotypes to be run in a job
# for testing 
n_chunk = 100

CHUNKS = range(1,n_chunk + 1)
fusion_dir = "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/"
ld_ref = fusion_dir + "LDREF"

prefix = out_folder + data_code


rule all:
    input:
        prefix + "/WEIGHTS/" + data_code + ".pos"

# take a phenotype bed file, split into a number of chunks
rule split_pheno:
    input: 
        pheno_df
    output:
        expand(prefix + "/pheno/" + data_code + ".pheno.{CHUNK}.bed", CHUNK = CHUNKS)
    params:
        out_prefix = prefix + "/pheno/" + data_code,
        script = "scripts/split_pheno.R"
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {input} -o {params.out_prefix} -n {n_chunk} -m {mode}"


# transpose covariates
rule transpose_cov:
    input:
        cov_df
    output:
        prefix + "/cov/" + data_code + ".covariates.tsv"
    params:
        script = "scripts/transpose_cov.R"
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {input} -o {output}"

# convert VCF to Plink and split by chrom
rule split_vcf:
    input:
        vcf
    output:
        expand( prefix + "/geno/" + data_code + ".{CHR}.bed", CHR = range(1,23) )
    run:
        for chromo in range(1,23):
            shell("ml plink; plink --vcf {input} --chr {chromo} --make-bed --out {prefix}/geno/{data_code}.{chromo} --extract {ld_ref}/1000G.EUR.{chromo}.bim")


# run Fahri/Gusev's bash script on each chunk
rule compute_twas_weights:
    input:
        pheno = prefix + "/pheno/" + data_code + ".pheno.{CHUNK}.bed",
        cov = prefix + "/cov/" + data_code + ".covariates.tsv",
        geno = expand( prefix + "/geno/" + data_code + ".{CHR}.bed", CHR = range(1,23) )
    output:
        prefix + "/logs/chunk.{CHUNK}.log"
    params:
        script = "scripts/compute_TWAS_functional_weights.sh",
        fusion_script = "scripts/efficient_FUSION.compute_weights.R",
        n1_chunk = n_chunk + 1,
        out_prefix = data_code,
        n_threads = 1
    shell:
        "sh {params.script} {input.pheno} {prefix}/geno/{data_code} {ld_ref} {input.cov} {params.n1_chunk} {params.out_prefix} {out_folder}/{data_code} {params.n_threads} {params.fusion_script} > {output}"

# create POS file listing the models
rule prepare_pos_file:
    input:
        expand( prefix + "/logs/chunk.{CHUNK}.log", CHUNK = CHUNKS)
    output:
        wgt_list = prefix + "/" + data_code + ".weights.list",
        wgt_profile = prefix + "/" + data_code + ".weights.profile",
        pos_file = prefix + "/WEIGHTS/" + data_code + ".pos"
    params:
        profile_script = fusion_dir + "/utils/FUSION.profile_wgt.R",
        prep_script = "scripts/create_pos_file.R"
    shell:
        "ml R/3.6.0;"
        "printf '%s\\n' {prefix}/WEIGHTS/*RDat > {output.wgt_list}; "
        "Rscript {params.profile_script} {output.wgt_list} > {output.wgt_profile}; "
        "Rscript {params.prep_script} --n_sample {n_sample} --weight_folder {prefix}/WEIGHTS --panel {data_code}"



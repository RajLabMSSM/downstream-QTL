import os

import pandas as pd


# for a QTL dataset, create FUSION weights for TWAS
# inputs (variables for each dataset)
# VCF - the full genotype matrix for common variants
# phenotype_df - this is produced by the QTL-mapping pipeline for expression or splicing phenotypes
# covariate_df - known covariates file - produced by QTL-mapping pipeline transposed to each covariate is a separate column
# data_code - the name of the dataset

# test data - the thoracic spinal cord eQTLs
data_code = "test"
vcf = "/sc/arion/projects/als-omics/WGS_QC/NYGC_Freeze02_European_Feb2020/WGS_QC_Pipeline/NYGC_Freeze02_European_Feb2020/output/chrAll_QCFinished_MAF0.01.anno.vcf.gz"
mode = "eQTL"
pheno_df = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/ThoracicSpinalCord_expression/ThoracicSpinalCord_expression.phenotype.tensorQTL.gene.bed.gz"
cov_df = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/ThoracicSpinalCord_expression/peer10/ThoracicSpinalCord_expression_peer10.gene.combined_covariates.txt"
out_folder = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/TWAS/TWAS_weights/"
# hardcoded inputs
# LD reference - from FUSION
# number of chunks - how many phenotypes to be run in a job
n_chunk = 400
CHUNKS = range(1,n_chunk + 1)
fusion_dir = "/sc/arion/projects/ad-omics/data/software/fusion_twas-master/"
ld_ref = fusion_dir + ""


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
        script = "scripts/split_pheno.R"
    shell:
        "ml R/3.6.0;"
        "Rscript {params.script} -i {input} -n {n_chunk} -m {mode}"


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
            shell("ml plink/1.9; plink --vcf {input} --chr {chromo} --make-bed --out {prefix}/geno/{data_code}.{chromo} ")


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
        n1_chunk = n_chunk + 1,
        out_prefix = data_code,
        n_threads = 1
    shell:
        "sh {params.script} {input.pheno} {ld_ref} {input.cov} {params.n1_chunk} {params.out_prefix} {out_folder} {params.n_threads} > {output}"

rule prepare_pos_file:
    input:
        expand( prefix + "/logs/chunk.{CHUNK}.log", CHUNK = CHUNKS)
    output:
        wgt_list = prefix + "/" + data_code + ".weights.list",
        wgt_profile = prefix + "/" + data_code + ".weights.profile",
        pos_file = prefix + "/WEIGHTS/" + data_code + ".pos"
    params:
        profile_script = fusion_dir + "/utils/FUSION.profile_wgt.R",
        prep_script = "scripts/preparePosFile.R"
    shell:
        "ml R/3.6.0;"
        "printf '%s\n' {prefix}/WEIGHTS/*RDat > {output.wgt_list}; "
        "Rscript {params.profile_script} {output.wgt_list} > {output.wgt_profile}; "
        "Rscript {params.prep_script} --phenotype {pheno_df} --prefix {data_code} --list {output.wgt_list} --output {output.pos_file}"

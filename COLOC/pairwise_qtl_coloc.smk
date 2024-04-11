import sys
import pandas as pd
# Pairwise QTL COLOC pipeline

outFolder = config["outFolder"]
inFolder = config["inFolder"]

target_files = [ x for x in glob.glob(inFolder + "/*genewide_QTL_pairs*tsv.gz" )]



shell.prefix("export PS1=""; ml anaconda3; CONDA_BASE=$(conda info --base); source $CONDA_BASE/etc/profile.d/conda.sh; ml purge; conda activate snakemake;")
rule all:
    input:
        #expand(outFolder + "{GWAS}/{QTL}/{QTL}_{GWAS}_coloc_results.snp.1kgp3_ld.tsv.gz", GWAS = gwas_data, QTL = qtl_data)
        #outFolder + "all_COLOC_results_merged_H4_0_no_LD.tsv.gz",
        #outFolder + "all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz",
        #outFolder + "pairwise_QTL_results.tsv.gz"
        #expand( outFolder + "{GWAS}/{QTL}_{GWAS}_COLOC.tsv", GWAS = gwas_data, QTL = qtl_data )
        #outFolder + "all_COLOC_summary_results.tsv.gz"

rule run_COLOC:
    output:
        outFolder + "{GWAS}/{QTL}_{GWAS}_COLOC.tsv",
        outFolder +  "{GWAS}/{QTL}_{GWAS}_COLOC.RData"
    params:
        script = "scripts/run_COLOC.R"
    shell:
        "mkdir -p {outFolder}{wildcards.GWAS};"
        "ml R/4.0.3;" #ml arrow;"
        "Rscript {params.script} -g {wildcards.GWAS} -q {wildcards.QTL} -o {outFolder}{wildcards.GWAS}/ "
   
rule merge_COLOC_snp_level:
    input:
        expand( outFolder + "{GWAS}/{QTL}_{GWAS}_COLOC.tsv", GWAS = gwas_data, QTL = qtl_data )
    output:
        outFolder + "all_COLOC_results_merged_H4_0_no_LD.tsv.gz",
        outFolder + "all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz"
    params:
        script = "scripts/merge_COLOC_results.R"
    shell:
        "ml R/3.6.0;" # LDLINKR only installed on 3.6
        "Rscript {params.script} --threshold 0 --inFolder {outFolder} --geneMeta {geneMeta};" 
        "Rscript {params.script} --threshold 0.5 --ld --inFolder {outFolder} --geneMeta {geneMeta} "

rule create_LD_tables:
    input:
        outFolder +  "{GWAS}/{QTL}_{GWAS}_COLOC.RData"
    output:
        outFolder +  "{GWAS}/{QTL}/{QTL}_{GWAS}_coloc_results.snp.1kgp3_ld.tsv.gz"
    params:
        script = "scripts/prepare_gwas_qtls.R",
        out = outFolder +  "{GWAS}/{QTL}/"
    shell:
       "conda activate echoR;"
       "Rscript {params.script} -g {wildcards.GWAS} -q {wildcards.QTL} --inFile {input} --outFolder {params.out}" 

rule prepare_QTL_COLOC:
    input:
        outFolder + "all_COLOC_results_merged_H4_0.5_with_LD.tsv.gz"
    output:
        outFolder + "pairwise_QTL_COLOC_targets.tsv.gz"
    params:
        script = "scripts/prepare_pairwise_qtl_coloc.R"
    shell:
        "ml R/4.2.0;"
        "Rscript {params.script} -i {input} -o {output} -m across -p 0.8"

rule run_QTL_COLOC:
    input:
        outFolder + "pairwise_QTL_COLOC_targets.tsv.gz"
    output:
        outFolder + "pairwise_QTL_results.tsv.gz"
    params:
        script = "scripts/run_COLOC.R"
    shell:
        "ml R/4.2.0;"
        "Rscript {params.script} -t {input} -o {output}"

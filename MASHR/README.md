## mashR
#### April 27, 2020

This repo is based on the GTEx mashR approach published by [Urbut et al, 2019](https://www.nature.com/articles/s41588-018-0268-8). It is also a good idea to take a look at the [eQTL outline page](https://stephenslab.github.io/mashr/articles/eQTL_outline.html) to get a sense of the pipeline. An important message is that mash needs to learn the sharing patterns from the data, so it is essential to use the nominal files from QTL analysis (NOT only the top results from the permuted file).

The goal of this page is to write step by step, how we can apply mashR in the RajLab datasets using the singularity system on Chimera server. I was able to run for the microglia cohort utilising the output files from QTL-mapping-pipeline. I will explicit the details here. 

First of all, the eQTL summary statistics must have the following columns:

gene_id 
variant_id      
tss_distance    
ma_samples      
ma_count        
maf     
pval_nominal    
slope   
slope_se

[Click here](https://github.com/stephenslab/gtexresults/blob/master/workflows/fastqtl_to_mash.ipynb) for more details. 

***************************************
Step by step to create the input files 
***************************************

1 - Convert parquet files to txt.gz. If you have only parquet files you can use the script [prepare_nominal_files.R](https://rajlabmssm.github.io/mashR/prepare_nominal_files.R) with 2 parameters: the path with the parquet files and the output filename. For example: 

```
ml R/3.6.2
./prepare_nominal_files.R ~/pd-omics/katia/Microglia/mic_314s/test_files THA_eur_expression_test.txt
```

2 - Create a list with ensembls to be used without a header. This step is fundamental to speed up the analysis. For the microglia project, we used ~19k genes. 

```
ENSG00000000003.14
ENSG00000000419.12
ENSG00000000457.14
ENSG00000000460.17
ENSG00000000938.13
``` 

3 - Create a file with the paths for the nominal files (named here as "QTLSumStats.list"). However, the Singularity will do NOT see the files outside his small world. So, after cloning the git, I copied these files to a folder inside the Singularity.

```
cat ~/pd-omics/katia/Microglia/mic_314s/input_mashR/QTLSumStats.list
STG_eur_expression_peer10.cis_qtl_nominal.txt.gz
THA_eur_expression_peer10.cis_qtl_nominal.txt.gz
SVZ_eur_expression_peer5.cis_qtl_nominal.txt.gz
MFG_eur_expression_peer10.cis_qtl_nominal.txt.gz
```

4 - Clone the Github page with GTEx results. It will create a folder that includes the scripts needed  to run the tool. 

```
cd ~/pd-omics_hydra/katia/softwares/ 
mkdir -p mash
cd mash
git clone https://github.com/stephenslab/gtexresults.git

cd gtexresults
ml singularity/3.5.2
singularity pull docker://gaow/hdf5tools
alias fastqtl2mash='singularity run -H $PWD hdf5tools_latest.sif'
```

5 - If you want to run an example you can try the code below. The results will be saved in the folder fastqtl_to_mash_output.

```
fastqtl2mash sos run workflows/fastqtl_to_mash.ipynb --data_list data/fastqtl/FastQTLSumStats.list --gene_list data/fastqtl/GTEx_genes.txt -j 8
```

6 - Now that we have the whole environment, we can start to substitute the files. I copied the nominal files, the QTLSumStats.list and the ensembl_list.txt for the folder mash/gtexresults/data/fastqtl. 

7 - Inside the main folder (mash/gtexresults/) you can run the sos workflow (I recommend submitting a job for that). 

```
cd ~/pd-omics_hydra/katia/softwares/mash/gtexresults/

bsub -q long -n 8 -W 144:00 -R "rusage[mem=20000]" -R "span[hosts=1]" -P acc_ad-omics -Ip /bin/bash

ml singularity/3.5.2
alias fastqtl2mash='singularity run -H $PWD hdf5tools_latest.sif'
fastqtl2mash sos run workflows/fastqtl_to_mash.ipynb --data_list data/fastqtl/QTLSumStats.list --gene_list data/fastqtl/ensembl_list.txt -j 8
```

Great! Now you have the "mashable" files to use as input. You will find a [description of the data here](https://stephenslab.github.io/gtexresults/gtexdata.html). 

Let's continue. Check the integrity of the data. Check if you get the same number of groups (genes) at the end of HDF5 data conversion: 

```
zcat THA_eur_expression_peer10.cis_qtl_nominal.txt.gz | cut -f1 | sort -u | wc -l

h5ls THA_eur_expression_peer10.cis_qtl_nominal.txt.h5 | wc -l
```

Now, you can run an implementation of the mashr eQTL workflow using an executable package. [Authors page is here](https://github.com/stephenslab/gtexresults/blob/master/workflows/mashr_flashr_workflow.ipynb).

***************************************
Step by step to run mashR
***************************************
To run mashR you will need a different docker image. So, you first pull the image and create an alias for that. 

```
cd /sc/hydra/projects/pd-omics/katia/softwares/mash/gtexresults 
ml singularity/3.5.2
singularity pull docker://gaow/mash-paper
alias mash='singularity run -H $PWD mash-paper_latest.sif'
```

Different pipelines can be run from the sos workflow. Each one is described below. 

Run the flash pipeline 
```
mash sos run workflows/mashr_flashr_workflow.ipynb flash --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
```

Set up the data-driven and canonical matrices 
```
mash sos run workflows/mashr_flashr_workflow.ipynb prior --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
```

Create the Vhat matrix of correlation 
```
mash sos run workflows/mashr_flashr_workflow.ipynb vhat_identity --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
mash sos run workflows/mashr_flashr_workflow.ipynb vhat_simple --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
mash sos run workflows/mashr_flashr_workflow.ipynb vhat_mle --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
```

Fit mash model (estimate mixture proportions)

Now we fit mash to the random tests using both data-driven and canonical covariances. NOTE from the author: "We have to fit using a random set of tests, and not a dataset that is enriched for strong tests".
```
mash sos run workflows/mashr_flashr_workflow.ipynb mash --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
```

Compute posterior summaries
```
mash sos run workflows/mashr_flashr_workflow.ipynb posterior --data fastqtl_to_mash_output/QTLSumStats.mash.rds -j 8
```

And now what? Now, you will be able to open the rds files on R and play! 
Example:
```
work_dir = "/sc/hydra/projects/pd-omics/katia/softwares/mash/gtexresults/mashr_flashr_workflow_output/"

mash_posterior = readRDS(paste0(work_dir, "QTLSumStats.mash.EZ.FL_PC3.V_simple.posterior.rds"))
names(mash_posterior)

```


> by Katia Lopes

>  RajLab





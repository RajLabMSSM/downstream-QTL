## mashR
#### July 17, 2019

> All code is in R and depends on packages from CRAN and/ or Bioconductor.

[Click here to see mashR for SV-eQTL](https://rajlabmssm.github.io/mashR/mashR.html)

Paper from [Urbut et al, 2019](https://www.nature.com/articles/s41588-018-0268-8)

[eQTL analysis outline](https://stephenslab.github.io/mashr/articles/eQTL_outline.html) by Mattew Stephens. 

Our 1st strategy was to analyze the restuls from sc-eQTL across 8 different cell types but, unfortunatelly it didn't work. We realized that we don't have enough overlap of gene-SNP pairs across conditions. For the GTEx analysis from Urbut et. al, the authors have used: "the strong data contained about 16k tests (the top eQTL per gene), and for the random data we used 20k randomly-selected tests. (If you suspect true effects are very sparse then you might want to increase the size of the random subset, say to 200k)." Therefore for this example, we've used a dataset of SV-eQTLs from 4 brain regions (Raj Lab unpublished). 

## Analysis outline

The basic analysis strategy for eQTL can be: 

1 - Check if your dataset have a good representative data overlapping across conditions. The mashR needs to learn from your data! 

2 - As input it's necessary to create two different matrices:

	* Bhat = In our case, we've used the slope data from fastQTL output. 
	* Shat = Standard error matrix. 

3 - mashR will learn correlation structure among null tests using random test.

4 - mashR will learn data-driven covariance matrices using strong tests.

5 - Fit the mashr model to the random tests. 

6 - Compute posterior summaries on the strong tests, using the model fit from step 4. 






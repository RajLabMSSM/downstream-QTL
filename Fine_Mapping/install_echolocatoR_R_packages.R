#options(echo=TRUE)

## INSTALLING R PACKAGES FOR ECHOLOCATOR
message(" * installing R packages required for echolocatoR")

# Regular R packages ------------------------------------

# current CRAN version of foreign (dependency for ChIPQC) needs R >= 4.0 - so specify legacy version

r_packages <- c("r.utils", "reticulate", "pbmcapply", "plotly","cowplot", "patchwork", "ggrepel", "curl", "gaston", "tidyverse", "BiocManager", "crayon", "roxygen2", "coloc", "haploR", "doBy")

required_packages <- r_packages[ ! r_packages %in% installed.packages() ]

for(lib in required_packages){
    install.packages(lib,dependencies=TRUE)
}


# several packages are no longer available on CRAN - get the last approved versions

if( ! "foreign" %in% installed.packages() ){
    install.packages("https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-76.tar.gz", dependencies = TRUE)
}

if( ! "XGR" %in% installed.packages() ){
    install.packages("https://cran.r-project.org/src/contrib/Archive/XGR/XGR_1.1.7.tar.gz", dependencies = TRUE)
}

if( ! "refGenome" %in% installed.packages() ){
    install.packages("https://cran.r-project.org/src/contrib/Archive/refGenome/refGenome_1.7.7.tar.gz", dependencies = TRUE)
}


# library(cowplot)
# library(ggrepel)
# library(curl)
# library(gaston)
# library(tidyr)
# library(BiocManager)
# library(refGenome)
# library(crayon)
# library(roxygen2) #roxygenize()


# library(coloc)
# install.packages("haploR", dependencies = TRUE)
# library(haploR)
# library(XGR)# install.packages("XGR")


# Bioconductor -------------------------------------

library(BiocManager)
bioc_packages <- c("supraHex", "graph", "Rgraphviz", "dnet", "GeneOverlap", "rtracklayer", "biomaRt", "Rsamtools", "snpStats", "bigsnpr", "ChIPQC")


# CHIPQC failing due to ERROR: dependency ‘foreign’ is not available for package ‘Hmisc’

required_packages <- bioc_packages[ ! bioc_packages %in% installed.packages() ]

for(lib in required_packages){
    BiocManager::install(lib)
}




# BiocManager::install(c("supraHex","graph","Rgraphviz","dnet"))
# library(GeneOverlap) #BiocManager::install("GeneOverlap")
# library(rtracklayer) #BiocManager::install("rtracklayer")
# library(biomaRt) # BiocManager::install("biomaRt")
# library(Rsamtools) # BiocManager::install("Rsamtools")
# library(snpStats) #BiocManager::install("snpStats")
# library(bigsnpr) # BiocManager::install("bigsnpr")
# library(ChIPQC); BiocManager::install('ChIPQC')



# github packages ------------------------------------
library(devtools)

install_github('jimhester/knitrBootstrap')
install_github("stephenslab/susieR")
install_github("variani/finemapr")
install_github("boxiangliu/locuscomparer")

# locusComparer requires a MySQL database - 

# library(fGWAS)
# *** SUSIE ****
# library(knitrBootstrap) #install_github('jimhester/knitrBootstrap')
# library(susieR) # devtools::install_github("stephenslab/susieR")

# *** finemapr ****
## finemapr contains: finemap, CAVIAR, and PAINTOR
# library(finemapr) # devtools::install_github("variani/finemapr")
# *** locuscomparer ****
# https://github.com/boxiangliu/locuscomparer
# library(locuscomparer); 

#devtools::install_github("boxiangliu/locuscomparer")

#install_fGWAS <- function(){
#   devtools::install_github("wzhy2000/fGWAS/pkg")
#   system("git clone https://github.com/wzhy2000/fGWAS.git; cd fGWAS; R CMD INSTALL pkg")
#}

#install_fGWAS()


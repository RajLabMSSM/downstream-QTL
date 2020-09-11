# Fine-Mapping

Scripts to accompany Brian Schilder's **echolocatoR** Fine-Mapping pipeline.

## Installing echolocatoR

first clone the echolocatoR repo

create a new conda environment using Brian's YAML file

```
conda env create -f echolocatoR/inst/conda/echoR.yml
```


activate the environment

run install_echolocator_R_packages.R to install all the R package dependencies

finally install echolocatoR while inside the environment: 

```
cd <folder where you cloned the repo, not inside the repo folder>
R CMD INSTALL echolocatoR/
```

## Fine-Mapping using COLOC pipeline results

to be described

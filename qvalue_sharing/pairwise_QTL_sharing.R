library("tidyverse")
library(purrr)


source_file <- "NYGC_sQTL_permutation_table.tsv"
target_file <- "NYGC_sQTL_nominal_table.tsv"

outFolder <- "NYGC_sQTLs/"

if( !dir.exists(outFolder) ){dir.create(outFolder) }

prepare_cmd <- function(source_path,target_path, outFolder, source_name, target_name){
    outFile <- paste0(outFolder, "/", source_name, "_", target_name, "_pvalues.tsv")
     cmd <- paste0("source activate QTL-pipeline; module purge; python prepare_qvalue.py -s ", source_path, " -t ", target_path, " -o ", outFile )
    print(cmd)
    if(file.exists(outFile) ){
        return("echo 0") }
    else{
        return(cmd)
    }
}

qvalue_cmd <- function(source_path, target_path, outFolder, source_name, target_name){
    inFile <- paste0(outFolder, "/",source_name, "_", target_name, "_pvalues.tsv")
    outFile <- paste0(outFolder, "/",source_name, "_", target_name, "_pi0.tsv")
    cmd <- paste0("source activate QTL-pipeline; Rscript run_qvalue.R -i ",  inFile, " -o ", outFile, " --source_name ", source_name, " --target_name ", target_name )
    print(cmd)
    if(file.exists(outFile) ){
        return("echo 0") }
    else{
        return(cmd)
    }
}

options(echo = TRUE)
# here each source_path is a permutation file
source_paths <- read_tsv(source_file)

# each target_path is a directory containing parquet files
target_paths <- read_tsv(target_file)

# for testing

#source_paths <- head(source_paths, 1 )
#target_paths <- head(target_paths,2)
#source_paths
#target_paths

for(i in 1:nrow(source_paths) ){
    # compare each target_path with each source_path
    for( j in 1:nrow(target_paths) ){
        system( prepare_cmd(source_path = source_paths$path[i], target_path = target_paths$path[j], outFolder = outFolder, source_name = source_paths$name[i], target_name = target_paths$name[j]  ) )
        system( qvalue_cmd(source_path = source_paths$path[i], target_path = target_paths$path[j], outFolder = outFolder, source_name = source_paths$name[i], target_name = target_paths$name[j] ) )
    }

}

results_files <- list.files(outFolder, pattern = "_pi0.tsv", full.names = TRUE)

all_res <- map_df(results_files, ~{read_tsv(.x) }) 

write_tsv(all_res, paste0(outFolder, "all_results.tsv" ) )

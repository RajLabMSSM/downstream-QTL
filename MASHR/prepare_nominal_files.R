#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
if (length(args)!=2) {
  stop("Usage: Rscript prepare_nominal_files.R parquet_input_folder output_file", call.=FALSE)
} 
input_folder = args[1]
output_file = args[2]

# install.packages("arrow")
library(arrow)

setwd(input_folder)

file_list <- list.files(path=input_folder, pattern = "*.parquet")

for (i in 1:length(file_list)){
  chr_spec <- read_parquet(file_list[i])
  file_name <- gsub("parquet", "txt", file_list[i])
  write.table(chr_spec, file = file_name, sep = "\t", row.names = F, quote = F)
}

system(paste0("gzip *.chr*.txt"))
system(paste0("zcat *.chr*.txt.gz > ", output_file))






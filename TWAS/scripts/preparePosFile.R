#!~/miniconda3/bin/Rscript
#Quickly prepare a pos file for TWAS, written by Fahri Kucukali
library(argparser, quietly=TRUE)
library(data.table)
library(dplyr)

p <- arg_parser("For creating a pos file for TWAS. Needs to have preparePosFile.sh in the path. Written by Fahri Kucukali. Version 1.00a")
p <- add_argument(p, "--phenotype", help="Molecular phenotype matrix file", nargs=1)
p <- add_argument(p, "--prefix", help="RDat file prefix", nargs=1)
p <- add_argument(p, "--list", help="List file", nargs=1)  
p <- add_argument(p, "--output", help="Output name for the pos file", nargs=1)  

script <- "scripts/preparePosFile.sh"

arg <- parse_args(p)

phenotype <- arg$phenotype

prefix <- arg$prefix

listname <- arg$list

output <- arg$output

system(paste0("bash ", script,phenotype, " " ,prefix))

guide <- fread(paste0(prefix,(".all_possible_WGT_GNAME_CHR_P0_P1.tsv")),header=FALSE)

listtable <- fread(listname, header=FALSE)

leftjoined <- left_join(listtable,guide,by="V1")

write.table(leftjoined,output,sep="\t",col.names=FALSE,quote=F,row.names=F)

system("echo -e 'WGT ID CHR P0 P1' > header_pos.txt")

system(paste0(("cat header_pos.txt "),output,(" > withheader."),output))

system(paste0(("cat withheader."),output,(" > "),output))

system(paste0(("rm withheader."),output))

system("rm header_pos.txt")

system(paste0(("rm "),prefix,(".all_possible_WGT_GNAME_CHR_P0_P1.tsv")))


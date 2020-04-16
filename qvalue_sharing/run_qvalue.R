library(qvalue)
library(optparse)
library(readr)

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file', default = ""),
        make_option(c('-s', '--sourceFile'), help= "the path to the source file" ),
        make_option(c('-t', '--targetDir'), help = "the path to the directory containing target files" ),
        make_option(c('--sourceName'), help = "the name of the source P-value distribution", default = "source"),
        make_option(c('--targetName'), help = 'the name of the target P-value distribution', default = "target"),
        make_option(c('--sourceType'), help = "whether the source file is eQTL or sQTL", default = "eQTL"),
        make_option(c('--targetType'), help = "whether the target files are eQTL or sQTL", defaul = "eQTL")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
sourceFile <- opt$sourceFile
targetFile <- opt$targetFile
sourceName <- opt$sourceName
targetName <- opt$targetName
sourceType <- opt$sourceType
targetType <- opt$targetType

library(dplyr)
library(arrow)

extract_source <- function(sourceFile, sourceType = "tensorQTL",mode ){
    source_df <- readr::read_tsv(sourceFile)
    if( sourceType == "tensorQTL"){
        varCol <- "variant_id"
        pCol <- "qval"
    }
    source_sig <- filter(source_df, qval < 0.05)
    # select gene, variant, p-value
    if(mode == "eQTL"){geneCol <- "phenotype_id"}
    if(mode == "sQTL"){geneCol <- "group_id"}
    source_sig <- select(source_sig, geneCol, varCol, pCol)
    names(source_sig) <- c("gene", "variant", "source") 
   return(source_sig)
}

extract_target_tensorQTL <- function(source_sig, targetDir, mode){
    targetFiles <- list.files(targetDir, ".parquet", full.names = TRUE)
    res <- purrr::map_df(targetFiles, ~{
        message("reading ", .x)
        target_df <- arrow::read_parquet(.x)
        pCol <- "pval_nominal"
        varCol <- "variant_id"
        if(mode == "eQTL"){geneCol <- "phenotype_id"}
        if(mode == "sQTL"){
            geneCol <- "group_id"
            target_df$group_id <- stringr::str_split_fixed(target_df$phenotype_id,":", n =5)[,5]
        }
        target_df <- select(target_df, gene = geneCol, variant = varCol, target = pCol)
        # get best variant per gene
        target_df_grouped <- group_by(target_df, gene) %>% summarise( target = min(target) )
        target_df_grouped <- left_join(target_df_grouped, target_df, by = c("gene", "target") )
        merged <- inner_join(source_sig, target_df, by = c("gene", "variant"))
        return(merged) 
    })

    return(res)
}


main <- function(){

message(" * running qvalue ")

#df <- read.table(inFile, header=TRUE, sep = "\t")

#stopifnot(c("source", "target") %in% names(df) )

# run qvalue


run_qvalue <- function(x){
    tryCatch(
        expr = {
            return(qvalue(x) )
            message("Successfully executed the log(x) call.")
        },
        error = function(e){
            message('Caught an error!')
            print(e)
        },
        warning = function(w){
            message('Caught an warning!')
            print(w)
        },
        finally = {
            message('All done, quitting.')
        }
    )    
}

res <- run_qvalue(df$target)

# if P-values are all highly significant then qvalue will throw an error
if("error" %in% class(res)){
    pi1 <- NA
}
if("pi0" %in% names(res) ){
    pi1 <- 1 - res$pi0
}

# create results table

res_df <- data.frame( source_name = sourceName, target_name = targetName, pi1 = pi1)

print(res_df)

write_tsv(res_df, path = outFile )  
}

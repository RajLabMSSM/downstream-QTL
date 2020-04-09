library(qvalue)
library(optparse)
library(readr)

option_list <- list(
        make_option(c('--outFile'), help='the path to the output file', default = ""),
         make_option(c('--sourceName'), help = "the name of the source P-value distribution", default = "source"),
        make_option(c('--targetName'), help = 'the name of the target P-value distribution', default = "target"),
    make_option(c('--inFile'), help = "the concatenated pvalue distribution table")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
inFile <- opt$inFile
sourceName <- opt$sourceName
targetName <- opt$targetName

df <- read.table(inFile, header=TRUE, sep = "\t")

stopifnot(c("source", "target") %in% names(df) )

# run qvalue

message(" * running qvalue ")

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


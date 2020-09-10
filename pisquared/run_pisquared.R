library(pisquared)

library(optparse)

option_list <- list(
        make_option(c('-o', '--outFolder'), help='the path to the output file', default = ".")
)

option.parser <- OptionParser(usage = "%prog [options] <pisquared input file>", description = "a script that combines QTL dataset into the format required by pisquared", option_list=option_list)
opt <- parse_args(option.parser, positional_arguments = 1)
input <- opt$args
options <- opt$options

outFolder <- options$outFolder

# read in pisquared input

df <- vroom::vroom(input)

df <- df[complete.cases(df),]

df <- df[ sample(1:nrow(df), size = 1e6, replace = FALSE),]

# use names to work out how many comparisons
comparisons <- names(df)[ 2:ncol(df) ]

all_comp <- expand.grid(comparisons, comparisons, stringsAsFactors = FALSE)

message(" * performing comparisons: ")
print(all_comp)

#all_comp <- all_comp[ all_comp$Var1 != all_comp$Var2, ]

# for testing
#all_comp <- all_comp[1,]

print(all_comp)

all_comp$pi2_jaccard  <- purrr::map_dbl(1:nrow(all_comp), ~{
    print(all_comp[.x,])
    p1 <- all_comp[.x, 1]
    p2 <- all_comp[.x, 2]
    pi2_res <- pi2_estimator(df[[p1]], df[[p2]], verbose = TRUE)
    return(pi2_res$jaccard)
})

outFile <- paste0(outFolder, "/pisquared_results.tsv")

readr::write_tsv(all_comp, path = outFile) 


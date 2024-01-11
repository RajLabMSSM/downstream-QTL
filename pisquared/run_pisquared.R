library(pisquared)

library(optparse)
set.seed(1234)
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

# randomly sample 1e6 rows
df <- df[ sample(1:nrow(df), size = 1e6, replace = FALSE),]

# use names to work out how many comparisons
comparisons <- names(df)[ 2:ncol(df) ]

all_comp <- expand.grid(comparisons, comparisons, stringsAsFactors = FALSE)

message(" * performing comparisons: ")
print(all_comp)

all_comp <- all_comp[ all_comp$Var1 != all_comp$Var2, ]


# for testing
#all_comp <- all_comp[1,]

print(all_comp)
# testing - run pisquared once and report all values
all_pisquared  <- purrr::map_df(1:nrow(all_comp), ~{
    print(all_comp[.x,])
    p1 <- all_comp[.x, 1]
    p2 <- all_comp[.x, 2]
    pi2_res <- pi2_estimator(df[[p1]], df[[p2]], verbose = TRUE, bin_pvalues=TRUE)
    res_list <- list(
        pi2_jaccard = pi2_res$jaccard,
        null_both = pi2_res$pi[1,1],       
        null1_alternate2 = pi2_res$pi[1,2],
        alternate1_null2 = pi2_res$pi[2,1],
        both_alternate = pi2_res$pi[2,2]
        )
    return(res_list)
} )

all_comp <- dplyr::bind_cols(all_comp, all_pisquared)

outFile <- paste0(outFolder, "/pisquared_results.tsv")

readr::write_tsv(all_comp, path = outFile) 


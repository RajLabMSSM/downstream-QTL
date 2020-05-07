# input:

# directory containing GWAS files
# BED file of GWAS hits to test
# directory containing QTL nominal results 

# process:
# take the GWAS hit coordinates
# expand by 1MB either side
# extract the GWAS P-value distribution using tabix

# read in QTL parquet file for that chromosome
# extract SNP-Gene pairs at the coordinates
# split pairs by Gene

# for each gene:
# get together the two distributions, merge by shared position
# run COLOC to get different colocalisations
# extract top QTL SNP for locus
# record results in table - GWAS locus, Gene, top QTL SNP, COLOC probabilities
library(purrr)
library(readr)
library(dplyr)
library(stringr)
library(arrow)
library(coloc)
library(optparse)

# flank coordinates by set number of bases (default 1MB)
# work on either a coordinate string or a dataframe containing chr start and end columns
makeCoords <- function(input, flank = 1e6){
    if(all(class(input) == "character")){
        coord <- parseCoords(input)
        coord$start <- coord$start - flank
        coord$end <- coord$end + flank
        coord <- paste0(coord$chr, ":", coord$start, "-", coord$end)    
        return(coord)
    }
    if( "data.frame" %in% class(input)){
        stopifnot( all(c("chr", "start", "end") %in% names(input) ) )
        stopifnot( flank >= 0)
        coords <- paste0(input$chr, ":", input$start - flank, "-",input$end + flank) 
        return(coords)
    }
}

parseCoords <- function(coords){
    split <- as.data.frame(stringr::str_split_fixed(coords, ":|-", 3), stringsAsFactors = FALSE)
    names(split) <- c("chr", "start","end")
    split$start <- as.numeric(split$start)
    split$end <- as.numeric(split$end)
    return(split)
}


# extract SNPs within coordinate range from a GWAS summary stat file
# account for different GWAS having different column naming and ordering with a GWAS_config.yaml file
# for each GWAS set the column numbers, the number of samples (N), the type of GWAS ("cc" or "quant") and the case proportion
# and whether the GWAS was hg19 or hg38
extractGWAS <- function(coord, config){
    # either read in config.yaml or Brian's CSV table
    config <- yaml::read_yaml(config) 
    
    gwas <- config$path
    chrCol <- config$chrCol
    posCol <- config$posCol
    pvalCol <- config$pvalCol
    betaCol <- config$betaCol
    seCol <- config$seCol
    # assume coord is a string following chr:start-end format
    chr <- parseCoords(coord)["chr"]
    
    cmd <- paste0( "ml bcftools; tabix ", gwas, coord ) 
    
    if( verbose == TRUE){
        print(cmd)
    }
    result <- as.data.frame(data.table::fread(cmd = cmd, nThread = 4) )
    #result <- read.table(text = system2(command = "tabix",args = cmd, stdout = TRUE, timeout = 1), sep = "\t", header = FALSE)
    message(" * GWAS extracted!")
    #result <- read.table(text = system(cmd, intern = TRUE), sep = "\t", header=FALSE)
    result$type <- GWAStype
    
    names(result)[chrCol] <- "chr"
    names(result)[posCol] <- "pos"
    names(result)[pvalCol] <- "pvalues"
    names(result)[betaCol] <- "beta"
    names(result)[seCol] <- "varbeta"
    # convert standard error to the variance
    result$varbeta <- result$varbeta^2
    names(result)[mafCol] <- "MAF"

    result$snp <-  paste0(result$chr, ":", result$pos)
    result <- dplyr::select( result, snp, pvalues, beta, varbeta, MAF)
    
    result <- as.list(result)
    
    result$N <- N
    result$type <- GWAStype
    
    if( GWAStype == "cc" ){
        result$s <- caseProp
    }
    return(result)
}

# lift over GWAS
# if GWAS mapped to hg19 then take GWAS SNPs and lift over to hg38
# column names are standardised already
liftOverGWAS <- function(gwas, chainPath = "/sc/hydra/projects/ad-omics/data/references/liftOver/hg19ToHg38.over.chain.gz"){
    # make genomicRanges object
    # liftOver using 
    gwas_gr <-  GenomicRanges::GRanges(
        seqnames=S4Vectors::Rle(gwas$chr),
      ranges = IRanges::IRanges(start=as.numeric(gwas$pos), end=as.numeric(gwas$pos)),
      strand = S4Vectors::Rle(rep('*',nrow(gwas)))
        )
    # lift over using rtracklayer
    # watch out for duplicate entries
    lifted_over <- rtracklayer::liftOver(gwas_gr, chainPath)
    stopifnot(length(lifted_over == length(gwas_gr) )
    gwas$chr <- seqnames(lifted_over)
    gwas$pos <- start(lifted_over)

    return(lifted_over)
}

# extract all nominal QTL P-values overlapping the flanked GWAS hit
# split by Gene being tested
# Nominal QTL associations are stored in parquet files, one for each chromosome
# random access isn't possible (yet) so you have to read in the entire file and subset out the region of interest
extractQTL <- function(coord, qtl = qtl_prefix, varCol = 2, geneCol = 1, pvalCol = 7, betaCol = 8, seCol = 9, mafCol =4, N = 100, sig_level = 0.05 ){
    coord_split <- parseCoords(coord)
    stopifnot(nrow(coord_split) == 1)
    chr <- unlist(coord_split["chr"])

    perm_file <- paste0(qtl, ".cis_qtl.txt.gz" )
    stopifnot(file.exists(perm_file) )
    
    # Read in permutation results
    # get out significant genes within locus coordinates
    perm_res <- readr::read_tsv(perm_file)
    perm_res <- dplyr::bind_cols(perm_res, parseCoords(perm_res$variant_id) )
    sig <- dplyr::filter(perm_res, qval < sig_level & chr ==  coord_split$chr & start >= coord_split$start & start <= coord_split$end )
    
    message( paste0( nrow(sig), " significant genes or splicing events at this locus" ) )
    if( nrow(sig) == 0 ){ return(NULL) }
    
    # read in nominal QTL associations
    parquet_file <- paste0(qtl, ".cis_qtl_pairs.", chr, ".parquet" )
    stopifnot(file.exists(parquet_file) )
    pq <- arrow::read_parquet(parquet_file)
    
    # default column numbers specify the relevant columns in the file - these may change over time 
    names(pq)[geneCol] <- "gene"
    
    pq <- dplyr::filter(pq, gene %in% sig$phenotype_id)
         
    variant_pos <- pq[[varCol]]
    pq$chr <- chr
    pq$pos <- unlist(parseCoords( variant_pos )["start"])
    pq$snp <- paste0(chr,':', pq$pos)
    names(pq)[pvalCol] <- "pvalues"
    names(pq)[betaCol] <- "beta"
    names(pq)[seCol] <- "varbeta"
    pq$varbeta <- pq$varbeta^2
    names(pq)[mafCol] <- "MAF"
    # retain only associations within locus coords
    pq_subset <- 
        dplyr::select(pq, gene, snp, pos, pvalues, beta, varbeta, MAF) %>%
        dplyr::filter( pos >= coord_split$start & pos <= coord_split$end)

    # split by gene, convert to list object
    pq_split <- 
        split(pq_subset, pq_subset$gene) %>%
        purrr::map( ~{ 
            x = as.list(.x)
            x$N <- N
            x$type <- "quant"
            return(x)
        })
    return(pq_split)
}

       
# running COLOC between a GWAS locus and all eQTLs within 1MB either side
# flank coordinates
# extract GWAS nominal SNPs from region
# find GWAS hit SNP for recording
# extract QTLs within range - at set significance level
# join together and run COLOC
# output a list of objects:
## COLOC results - this should include the GWAS locus, the gene of interest, the top QTL variant and the top QTL p-value
## object - this should be a table combining the two input datasets, inner joined on "snp", to be used for plotting.
runCOLOC <- function(gwas_prefix, qtl_prefix, hit){
    coord <- parseCoords(hit)
    #coord <- makeCoords(hit, flank = 0)
    range <- makeCoords(coord, flank = 1e6)
    # extract GWAS and QTL summary for given coord
    message(" * extracting GWAS")
    g <- extractGWAS(range, gwas = gwas_prefix)
    
    # if GWAS is hg19 then lift over to hg38
    if(GWAS.build == "hg19"){
        message(" * GWAS is hg19. Lifting to hg38")
        g <- liftOverGWAS(g)
    }

    # get hit info from GWAS
    hit_snp <- paste(parseCoords(hit)[1:2],collapse = ":")
    hit_info <- as.data.frame(g, stringsAsFactors = FALSE) %>% filter(snp == hit_snp) %>% select(GWAS_SNP = snp, GWAS_P = pvalues, GWAS_beta = beta, GWAS_MAF = MAF) 
    message(" * extracting QTLs")
    q <- extractQTL(range, qtl = qtl_prefix)
    if( is.null(q) ){ return(NULL) }
    # for each gene extract top QTL SNP    
    qtl_info <- q %>% 
        map_df( ~{ 
            as.data.frame(.x, stringsAsFactors=FALSE) %>% 
            arrange(pvalues) %>% head(1) %>% 
            select(gene, QTL_SNP = snp, QTL_P = pvalues, QTL_Beta = beta, QTL_MAF = MAF)
        }) 
    message(" * running COLOC")
    coloc_res <- 
        purrr::map_df( q, ~{
            as.data.frame(t(coloc.abf(dataset1 = g, dataset2 = .x)$summary), stringsAsFactors = FALSE)
        }, .id = "gene") %>% 
        dplyr::mutate(GWAS_SNP = paste(parseCoords(hit)[1:2],collapse = ":")  ) %>% 
        dplyr::select(gene, GWAS_SNP, everything() )
    # bind coloc_res to other tables
    full_res <- full_join(qtl_info,
            full_join(hit_info, coloc_res, by = "GWAS_SNP"), 
                by = "gene") %>% 
        select( gene, starts_with("GWAS"), starts_with("QTL"), everything() )
    return(full_res) 
}

option_list <- list(
        make_option(c('-o', '--outFile'), help='the path to the output file', default = ""),
        make_option(c('--hits'), help= "the path to a file containing GWAS hits" ),
        make_option(c('-p', '--gwas_prefix'), help = "the prefix of the GWAS files"),
        make_option(c('-q','--qtl_prefix'), help = "the directory containing QTL nominal results with the prefix of the dataset; eg /results/Brain_expression/peer30/Brain_expression_peer30", default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFile <- opt$outFile
hits_file <- opt$hits
qtl_prefix <- opt$qtl_prefix
gwas_prefix <- opt$gwas_prefix


#gwas_prefix <- "/sc/arion/projects/als-omics/ALS_GWAS/Nicolas_2018/processed/Nicolas_2018_processed_"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_expression/peer30/LumbarSpinalCord_expression_peer30"
#qtl_prefix <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/LumbarSpinalCord_splicing/peer20/LumbarSpinalCord_splicing_peer20"

#hits_file <- "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/downstream-QTL/COLOC/Nicolas_2018/Nicolas_2018_hits_1e-7.tsv"
#outFile <- "test_coloc_results.tsv"

options(echo = TRUE)


hits <- read_tsv(hits_file)
hit_coords <- makeCoords(hits, flank = 0)
# for testing
#hit_coords <- hit_coords[1]

all_res <- purrr::map_df(1:length(hit_coords), ~{
    message(hit_coords[.x])
    res <- runCOLOC(gwas_prefix, qtl_prefix, hit = hit_coords[.x])
    if(is.null(res) ){return(NULL) }
    return(res)
})


readr::write_tsv(all_res, path = outFile)

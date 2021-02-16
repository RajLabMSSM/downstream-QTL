# WRITTEN BY FAHRI KUCUKALI?
# TWEAKED BY JACK HUMPHREY

ml R/3.6.0
ml plink/1.9
#!/bin/sh
# MAKE SURE efficient_FUSION.compute_weights.R IS IN YOUR PATH
# FILL IN THESE PATHS
GCTA="/home/fahri/AMP_AD/TWAS/fusion_twas-master/gcta_nr_robust"
PLINK="/home/fahri/miniconda3/bin/plink"
GEMMA="/home/fahri/miniconda3/bin/gemma"
# ALTERNATIVELY: ENSURE THAT plink, gcta, gemma CAN BE CALLED FROM PATH AND REMOVE --PATH_* FLAGS BELOW
# PATH TO DIRECTORY CONTAINING LDREF DATA (FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)
LDREF=$3
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY

# PATH TO GEUVADIS GENE EXPRESSION MATRIX:
PRE_GEXP=$1
# GEUVADIS DATA WAS DOWNLOADED FROM https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/

# PATH TO PREFIX FOR GEUVADIS GENOTYPES SPLIT BY CHROMOSOME
# SUBSAMPLE THESE TO THE LDREF SNPS FOR EFFICIENCY
PRE_GENO=$2

# PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)

OUT_DIR_FAHRI=$7

OUT_DIR="$7/$PRE_GEXP.WEIGHTS"

# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir $OUT_DIR

# ROWS IN THE MATRIX TO ANALYZE (FOR BATCHED RUNS)
BATCH_START=1
BATCH_END=$5

# OUT NAME and THREAD ADDED by Fahri

OUTNAME=$6

THREAD=$8

# IDs

cat $PRE_GEXP | head -n1 | tr '\t' '\n' | tail -n+5 > $PRE_GENO.IDs

# COVARIATEs

COVAR=$4

# --- BEGIN SCRIPT:

NR="${BATCH_START}_${BATCH_END}"
mkdir --parents $PRE_GEXP.tmp/$NR
mkdir --parents $PRE_GEXP.hsq/$NR
mkdir --parents $PRE_GEXP.out/$NR

# Loop through each gene expression phenotype in the batch
cat $PRE_GEXP | awk -vs=$BATCH_START -ve=$BATCH_END 'NR > s && NR <= e' |  while read PARAM; do

# Get the gene positions +/- 1Mb
    CHR=`echo $PARAM | awk '{ print $1 }'`
    P0=`echo $PARAM | awk '{if($2 - 1e6 <0) print 1; else print $2 - 1e6}'`
    P1=`echo $PARAM | awk '{ print $3 + 1e6 }'`
    GNAME=`echo $PARAM | awk '{ print $4 }'`

    OUT="$PRE_GEXP.tmp/$NR/$PRE_GEXP.$GNAME"

    echo $GNAME $CHR $P0 $P1

    # Pull out the current gene expression phenotype
    echo $PARAM | tr ' ' '\n' | tail -n+5  > $OUT.intermediatepheno

    paste $PRE_GENO.IDs $PRE_GENO.IDs $OUT.intermediatepheno > $OUT.pheno

    # Get the locus genotypes for all samples and set current gene expression as the phenotype
    $PLINK --bfile $PRE_GENO.$CHR --allow-no-sex --pheno $OUT.pheno --make-bed --out $OUT --keep $OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1 --extract $LDREF/hwe1e6.1000G.EURn404.GRCh38_fk.chr$CHR.bim

    # Process all samples together (for reference purposes only since this is multi-ethnic data)
    FINAL_OUT="$OUT_DIR/$OUTNAME.$GNAME"

    Rscript -e "Sys.getenv('OPT_NUM_THREADS'); opt <- list(threads = 1); Sys.setenv('OPT_NUM_THREADS' = opt$threads); Sys.getenv('OPT_NUM_THREADS')"

    Rscript efficient_FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 2 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet --covar $COVAR --hsq_p 0.05 --threads $THREAD 

    Rscript -e "Sys.getenv('OPT_NUM_THREADS')"

# ALTERNATIVELY ADD COVARIATES HERE USING THE --covar FLAG
# MINIMAL COMMAND IS: `Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT`

# Append heritability output to hsq file
    cat $FINAL_OUT.hsq >> $PRE_GEXP.hsq/$NR.hsq

# Clean-up just in case
    rm -f $FINAL_OUT.hsq $OUT.tmp.*

    # Remove all intermediate files
    rm $OUT.*

# GO TO THE NEXT GENE
done

echo " * Run finished!"

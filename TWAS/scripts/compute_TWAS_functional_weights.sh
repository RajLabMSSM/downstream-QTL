# WRITTEN BY FAHRI KUCUKALI?
# TWEAKED BY JACK HUMPHREY

ml R/3.6.0
ml plink
ml gemma
#!/bin/sh
# MAKE SURE efficient_FUSION.compute_weights.R IS IN YOUR PATH
# FILL IN THESE PATHS
GCTA="/sc/arion/projects/ad-omics/data/software/fusion_twas-master/gcta_nr_robust"
PLINK="plink"
GEMMA="gemma"
#PLINK="/home/fahri/miniconda3/bin/plink"
#GEMMA="/home/fahri/miniconda3/bin/gemma"
# ALTERNATIVELY: ENSURE THAT plink, gcta, gemma CAN BE CALLED FROM PATH AND REMOVE --PATH_* FLAGS BELOW
# PATH TO DIRECTORY CONTAINING LDREF DATA (FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)
LDREF=$3
echo " * using LD reference files in $3"
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY

# PATH TO GEUVADIS GENE EXPRESSION MATRIX:
PRE_GEXP=$1
echo " * using phenotype matrix $1 "
# GEUVADIS DATA WAS DOWNLOADED FROM https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/

# PATH TO PREFIX FOR GEUVADIS GENOTYPES SPLIT BY CHROMOSOME
# SUBSAMPLE THESE TO THE LDREF SNPS FOR EFFICIENCY
PRE_GENO=$2
echo " * genotypes: $2 "
# PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)

OUT_DIR_FAHRI=$7
echo " * out directory: $7"

OUT_DIR="$7/WEIGHTS"

# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir -p $OUT_DIR

# ROWS IN THE MATRIX TO ANALYZE (FOR BATCHED RUNS)
BATCH_START=1
#BATCH_END=$5
BATCH_END=$(wc -l $PRE_GEXP | awk '{print $1}')

# OUT NAME and THREAD ADDED by Fahri

OUTNAME=$6
echo " * out file name: $6"
THREAD=$8

FUSION_SCRIPT=$9
echo " * using FUSION script to compute weights: $9"
# IDs

cat $PRE_GEXP | head -n1 | tr '\t' '\n' | tail -n+5 > ${PRE_GENO}.IDs

# COVARIATEs

COVAR=$4
echo " * covariates: $4"
# --- BEGIN SCRIPT:

NR="${BATCH_START}_${BATCH_END}"
mkdir --parents $PRE_GEXP.tmp/$NR
mkdir --parents $PRE_GEXP.hsq/$NR
mkdir --parents $PRE_GEXP.out/$NR

if [ -e ${PRE_GEXP}.hsq/${NR}.hsq ];then
    rm -f ${PRE_GEXP}.hsq/${NR}.hsq
fi

# Loop through each gene expression phenotype in the batch
cat $PRE_GEXP | awk -vs=$BATCH_START -ve=$BATCH_END 'NR > s && NR <= e' |  while read PARAM; do

PRE_GEXP_BASE=$(basename $PRE_GEXP)

# Get the gene positions +/- 1Mb
    CHR=`echo $PARAM | awk '{ print $1 }'`
    P0=`echo $PARAM | awk '{if($2 - 1e6 <0) print 1; else print $2 - 1e6}'`
    P1=`echo $PARAM | awk '{ print $3 + 1e6 }'`
    GNAME=`echo $PARAM | awk '{ print $4 }'`
    echo " * creating TWAS weight for $GNAME !"
    OUT="$PRE_GEXP.tmp/$NR/$PRE_GEXP_BASE.$GNAME"

    echo $GNAME $CHR $P0 $P1

    # Pull out the current gene expression phenotype
    echo $PARAM | tr ' ' '\n' | tail -n+5  > $OUT.intermediatepheno

    paste ${PRE_GENO}.IDs ${PRE_GENO}.IDs $OUT.intermediatepheno > $OUT.pheno

    # Get the locus genotypes for all samples and set current gene expression as the phenotype
    # extract with only 1000G SNPs is done already when VCF is converted to per-chrom plink files
    $PLINK --bfile $PRE_GENO.$CHR --allow-no-sex --pheno $OUT.pheno --make-bed --out $OUT --keep $OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1 #--extract $LDREF/1000G.EUR.$CHR.bim

    # Process all samples together (for reference purposes only since this is multi-ethnic data)
    FINAL_OUT="${OUT_DIR}/${OUTNAME}.${GNAME}"
    echo " * writing to $FINAL_OUT "
    Rscript -e "Sys.getenv('OPT_NUM_THREADS'); opt <- list(threads = 1); Sys.setenv('OPT_NUM_THREADS' = opt$threads); Sys.getenv('OPT_NUM_THREADS')"

    Rscript ${FUSION_SCRIPT} --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 2 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet --covar $COVAR --hsq_p 0.05 --threads $THREAD 

    Rscript -e "Sys.getenv('OPT_NUM_THREADS')"

# ALTERNATIVELY ADD COVARIATES HERE USING THE --covar FLAG
# MINIMAL COMMAND IS: `Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT`

# Append heritability output to hsq file
    cat ${FINAL_OUT}.hsq >> ${PRE_GEXP}.hsq/${NR}.hsq

# Clean-up just in case
    rm -f ${FINAL_OUT}.hsq ${OUT}.tmp.*

    # Remove all intermediate files
    rm $OUT.*

# GO TO THE NEXT GENE
done

echo " * Run finished!"

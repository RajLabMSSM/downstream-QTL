#!/bin/sh
PRE_GEXP=$1

RDAT_PREFIX=$2

# Loop through each gene expression phenotype in the batch
zcat $PRE_GEXP | cut -f1,2,3,4 | awk '{print $4":"$1":"$2":"$3}' | grep -v gene_id | while read PARAM; do

    # Get the gene positions +/- 1Mb
    CHR=`echo $PARAM | cut -f2 -d":"`
    P0=`echo $PARAM | cut -f3 -d":" | awk '{if($1 - 1e6 <0) print 1; else print $1 - 1e6}'`
    P1=`echo $PARAM | cut -f4 -d":" | awk '{ print $1 + 1e6 }'`
    GNAME=`echo $PARAM | cut -f1 -d":"`

echo ""$RDAT_PREFIX"."$GNAME".wgt.RDat" $GNAME $CHR $P0 $P1

done > $RDAT_PREFIX.all_possible_WGT_GNAME_CHR_P0_P1.tsv

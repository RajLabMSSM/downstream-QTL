
JARFILE=/sc/hydra/projects/ad-omics/data/software/METASOFT/Metasoft.jar
PTABLE=/sc/hydra/projects/ad-omics/data/software/METASOFT/HanEskinPvalueTable.txt
inFile=$1
outFile=$2
java -jar $JARFILE -pvalue_table $PTABLE -input $inFile -mvalue -mvalue_p_thres 1E-4 -output $outFile -verbose

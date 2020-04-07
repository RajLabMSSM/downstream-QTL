ml bcftools/1.9

GWAS=$1

# sort GWAS by coordinate and position
zless $GWAS | awk 'NR >1' | sort -k 3,3 -k4,4n - | bgzip > ${GWAS}_sorted.bgz


tabix -b 4 -s 3 -e 4 -f ${GWAS}_sorted.bgz


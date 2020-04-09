
# prepare q-value sharing

import pyarrow.parquet as pq
import pandas as pd
import glob
from optparse import OptionParser


# set arguments
# --source - source file of P-values
# --source_type - string specifying tool used to generate
# --target - source file or directory containing files
# --target_type - default 'parquet'

parser = OptionParser()

parser.add_option("-s", "--source", dest="source")
parser.add_option("-t", "--target", dest="target")
parser.add_option("--source_type", dest="source_type", default = "tensorQTL")
parser.add_option("--target_type", dest="target_type", default = "tensorQTL")
parser.add_option("-o", "--output_file", dest="outFile")
(options, args) = parser.parse_args()

source = options.source
target = options.target
source_type = options.source_type
target_type = options.target_type
outFile = options.outFile
# inputs:
## A permuted P-value list

#source = "Hippocampus_expression_peer15.cis_qtl.txt.gz"
#source_type = "tensorQTL"

#target = "."
#target_type = "tensorQTL"

#outFile = "test_qvalue_input.tsv"

print(" * preparing data for qvalue comparison ")

print(" * source dataset: " + source )

if( target_type == "tensorQTl" ):
    print(" * target files: " + target)
else:
    print(" * target file: " + target)

# functions
def readSource(source, source_type):
    ## read in permutation file 
    source_df = pd.read_csv(source, sep = "\t")
    
    if source_type == "tensorQTL":
        geneCol = 'phenotype_id'
        varCol = 'variant_id'
        pCol = 'qval'
    
    # test if the supposed columns are in the file

    # subset source df according to type of input table
    source_df = source_df[[geneCol, varCol, pCol]]

    source_df = source_df.rename(columns = {pCol : "source" })

    # threshold for significance
    source_sig = source_df[ source_df['source'] < 0.05 ]

    return(source_sig)

# if tensorQTL target to be used, find all parquet files in a directory
# read in and merge with source_sig
def target_tensorQTL(target, source_sig):
    ## read in nominal files
    #nom_dir = "/sc/arion/projects/als-omics/QTL/NYGC_Freeze02_European_Feb2020/QTL-mapping-pipeline/results/Cerebellum_expression/peer30"

    nom_files = glob.glob(target + "/*parquet")

    nominal_tables = []

    for nom in nom_files:
        print(" * reading : " + nom )
        target_df = pq.read_table(nom, columns = ['phenotype_id', 'variant_id', 'pval_nominal']).to_pandas()
        target_df = target_df.rename(columns={"pval_nominal": "target"})
        ## match the two together on shared SNP and Gene columns
        merge = pd.merge(source_sig,target_df, on = ["phenotype_id", "variant_id"], how = "inner") 
        nominal_tables.append(merge)


    # concatenate together
    print(" * merging files together")
    merged_df = pd.concat(nominal_tables)

    return(merged_df)

source_sig = readSource(source, source_type)

merged_df = target_tensorQTL(target,source_sig)

merged_df.to_csv(outFile, sep = "\t", header = True, index = False  )


import pandas as pd 
import sys
import ast
import os 
if len(sys.argv)!=3:
    sys.exit("USAGE: python run_benchmark_vs_decoys.py inputfile outputdir")

folderprefix="../clusters/data"
filesuffix="finalList.tsv"
inputfile=sys.argv[1]
outputdir=sys.argv[2]

data=pd.read_csv(inputfile,sep=";")

for index, row in data.iterrows():
    if index==0:
        accession=row["Accession"]

        decoy_candidates=ast.literal_eval(row["DecoyCandidates"])

        decoystring=str([f"{folderprefix}/{decoy}{filesuffix}" for decoy in decoy_candidates]).replace(" ","").replace("'","")[1:-1]
        # if not os.path.exists(f"{outputdir}/{accession}"):
        os.makedirs(f"{outputdir}/{accession}",exist_ok=True)    
        print(f"python benchmark_vs_decoys.py {folderprefix}/{accession}{filesuffix} {decoystring} 0")
        # print(f"python analyzeBenchmark2Enrichment.py output/{accession}/{accession}finalList_out.csv output/{accession}/{accession}")
    
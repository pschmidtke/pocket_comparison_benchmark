import pandas as pd 
import sys
import ast
import os 
if len(sys.argv)!=3:
    sys.exit("USAGE: python run_benchmark_vs_decoys.py inputfile outputdir")

folderprefix="../clusters/data"
filesuffix="clusterRepresentatives.csv"
inputfile=sys.argv[1]
outputdir=sys.argv[2]

data=pd.read_csv(inputfile,sep=";")

for index, row in data.iterrows():
    accession=row["Accession"]

    decoy_candidates=ast.literal_eval(row["DecoyCandidates"])

    decoystring=str([f"{folderprefix}/{decoy}{filesuffix}" for decoy in decoy_candidates]).replace(" ","").replace("'","")[1:-1]
    if not os.path.exists(f"{outputdir}/{accession}/{accession}_f1_score.png"):
        # print(f"{outputdir}/{accession}/{accession}_f1_score.png")
        print(f"python benchmark_vs_decoys.py {folderprefix}/{accession}{filesuffix} {decoystring} outputClusterRepresentatives/{accession}")
        print(f"rm -rf tmp/{accession}clusterRepresentatives")
        print(f"python analyzeBenchmark2Enrichment.py outputClusterRepresentatives/{accession}/{accession}clusterRepresentatives_out.csv outputClusterRepresentatives/{accession}/{accession}")
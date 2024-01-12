import sys
import pandas as pd
import ast 
import gemmi
import subprocess
import numpy as np
import os
if(len(sys.argv) != 4):
    print("Usage: python3 testmol2stuff.py actives.csv decoys1.csv,decoys2.csv,decoys3.csv outputdir")
    sys.exit()



# Replace input.cif with your file and specify the chain ID and output file name

fname=sys.argv[1]
decoyfnames=sys.argv[2].split(",")
outputdir=sys.argv[3]
decoys=[pd.read_csv(decoyfname, sep="\t") for decoyfname in decoyfnames if os.path.exists(decoyfname)]
decoysdf=pd.concat(decoys)
decoysdf["decoy"] = 1
basename=os.path.basename(fname).split(".")[0]
outname=f"{outputdir}/{basename}_out.csv"
df = pd.read_csv(fname, sep="\t")
df["decoy"] = 0
alldata=pd.concat([df,decoysdf])
cwd = os.getcwd()

pdblist=[]
if not os.path.exists("tmp/"+basename):
    os.makedirs("tmp/"+basename+"/pdb",exist_ok=True)
for index, row in alldata.iterrows():
    pdbid=row["pdbid"].lower()

    if not os.path.exists("tmp/"+basename+"/pdb/"+pdbid+".pdb"):
        residues = ast.literal_eval(row["residues"])
        # print(residues)
        if (residues[0] != None ):  #and pdbid=="1UY8"
            filepath=f"/mnt/share/data2/pdb/structures_cifs/{pdbid[1:3]}/{pdbid}.cif.gz"
            st = gemmi.read_structure(filepath)
            for chain in st[0]:
                for residue in chain:
                    if str(residue.subchain)+":"+str(residue.label_seq) in residues:
                        residue.flag="s"

            selection = gemmi.Selection().set_residue_flags('s')
            pocket_st=selection.copy_structure_selection(st)
            try:
                pocket_st.write_minimal_pdb(f"tmp/{basename}/pdb/{pdbid}.pdb")
            except Exception as e:
                print(f"Error writing {pdbid} to pdb: {e}")
                continue

subprocess.run(["/usr/bin/bash","Step0-cabbage.sh", f"tmp/{basename}/pdb/", f"tmp/{basename}/pdb/"],env={"PATH": os.environ["PATH"]+":/home/ubuntu/3decision/pc/pocketmatch/PocketMatch_v2.1/cabbage-file_maker"})
os.chdir(f"tmp/{basename}/pdb/")
subprocess.run(["Step3-PM_typeA","outfile.cabbage"],env={"PATH": os.environ["PATH"]+":/home/ubuntu/3decision/pc/pocketmatch/PocketMatch_v2.1"})
decoydict={}
os.makedirs(f"{cwd}/outputClusterRepresentatives/{basename.strip('clusterRepresentatives')}",exist_ok=True)
with open("PocketMatch_score.txt","r") as f:
    print(f"opening {cwd}/{outname}")
    with open(f"{cwd}/{outname}","w") as fout:
        fout.write("pdb_id1 pdb_id2 sim pdb1decoy pdb2decoy\n")
        lines=f.readlines()
        for lidx,line in enumerate(lines):
                tmp=line.split(" ")
                try:
                    pdbid1=tmp[0].split(".")[0]
                    pdbid2=tmp[1].split(".")[0]
                
                    score=tmp[3].strip()
                    if score=="NULL____":
                        score="0.0"
                    if pdbid1 not in decoydict: 
                        decoydict[pdbid1]=str(alldata.loc[alldata["pdbid"]==pdbid1.upper(),"decoy"].values[0])
                    if pdbid2 not in decoydict: 
                        decoydict[pdbid2]=str(alldata.loc[alldata["pdbid"]==pdbid2.upper(),"decoy"].values[0])
                    decoy1=decoydict[pdbid1]
                    decoy2=decoydict[pdbid2]
                    lines[lidx]=f"{pdbid1} {pdbid2} {score} {decoy1} {decoy2}\n"
                    fout.write(f"{pdbid1} {pdbid2} {score} {decoy1} {decoy2}\n")
                except Exception as e:
                    pass
                
                # fout.flush

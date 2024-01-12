import MDAnalysis as mda
import sys
import pandas as pd
import ast 
import gemmi
import subprocess
import numpy as np
import os
print(len(sys.argv))
if(len(sys.argv) != 4):
    print("Usage: python3 benchmark_vs_decoys.py actives.csv decoys1.csv,decoys2.csv,decoys3.csv blacklist")
    sys.exit()



# Replace input.cif with your file and specify the chain ID and output file name

fname=sys.argv[1]
decoyfnames=sys.argv[2].split(",")
blacklisted=sys.argv[3].split(",")
decoys=[pd.read_csv(decoyfname, sep="\t") for decoyfname in decoyfnames if os.path.exists(decoyfname)]
decoysdf=pd.concat(decoys)
decoysdf["decoy"] = 1
basename=os.path.basename(fname).split(".")[0].strip("finalList")
outname=f"output/{basename}/{basename}_out.csv"
blacklistname=f"output/{basename}/blacklisted.txt"
df = pd.read_csv(fname, sep="\t")
df["decoy"] = 0
alldata=pd.concat([df,decoysdf])

pdblist=[]
for index, row in alldata.iterrows():
    pdbid=row["pdbid"].lower()
    if not os.path.exists("tmp/mol2/"+pdbid+".mol2"):
        residues = ast.literal_eval(row["residues"])
        if (residues[0] != None ):  #and pdbid=="1UY8"
            # print(pdbid)
            chain_names = np.unique([ residue.split(":")[0] for residue in residues ])
            # print(chain_names)
            resnumbers= [ residue.split(":")[1] for residue in residues ]
            if (len(chain_names)==1):
                filepath=f"/mnt/share/data2/pdb/structures_cifs/{pdbid[1:3]}/{pdbid}.cif.gz"
                # try:
                st = gemmi.read_structure(filepath)
                sel = gemmi.Selection(f"//{chain_names[0]}/*/*")
                try:
                    ca_st = sel.copy_structure_selection(st)
                    for chain in ca_st[0]:
                        for residue in chain:
                            residue.seqid.num = residue.label_seq
                            residue.subchain = chain.name
                
                    ca_st.write_minimal_pdb(f"tmp/pdb/{pdbid}.pdb")
                except Exception as e:
                    print(f"Error writing {pdbid} to pdb: {e}")
                    continue
                subprocess.run(["obabel", f"tmp/pdb/{pdbid}.pdb", "-O", f"tmp/mol2/{pdbid}.mol2","--title",pdbid])
                with open(f"tmp/mol2/{pdbid}.mol2", "r") as f:
                    lines = f.readlines()
                    selected=[]
                    bonds=0
                    skip=0
                    skippedLines=[]
                    for lidx,line in enumerate(lines):
                        if(len(line)>60):
                            skip=0
                            # print("    375 CZ          0.5440  -18.2220    18.5150 C.cat    25 ARG85       0.0000 DICT")
                            # print(f"{line[0:8].strip():>7} {line[9:13]} {line[16:26].strip():>13}  {line[28:36].strip():>8}   {line[37:46].strip():>8} {line[47:53]}  {line[53:56]} {line[58:65]}{line[68:].strip():>11}")
                            selectedstring=f" 0.0000"
                            if(line[47:53].strip() == "H"):
                                skippedLines.append(lidx)
                            if(line[53:56].strip()) in resnumbers and line[7:12].strip() == "CA":
                                selectedstring="40.0000"
                            lines[lidx]=f"{line[0:8].strip():>7} {line[9:13]} {line[16:26].strip():>13}  {line[28:36].strip():>8}   {line[37:46].strip():>8} {line[47:53]}  {line[53:56]} {line[58:65]}{selectedstring:>11}\n"
                        elif(line.strip()=="@<TRIPOS>UNITY_ATOM_ATTR"):
                            skip=1
                        elif(line.strip()=="@<TRIPOS>BOND"):
                            bonds=1
                            skip=0
                        elif(bonds==1):
                            lines[lidx]=f"{line[0:7].strip():>6}{line[7:12].strip():>5}{line[12:18].strip():>5} {line[19:].strip():<4}\n"
                        if(skip==1):
                            skippedLines.append(lidx)
                    finalLines=[lines for x,lines in enumerate(lines) if x not in skippedLines]
                    # print(f"{line[0:7].strip():>6}{line[7:12].strip():>5}{line[12:18].strip():>5} {line[19:].strip():<4}")
                with open(f"tmp/mol2/{pdbid}.mol2", "w") as f:
                    f.writelines(finalLines)
                    pdblist.append(pdbid)
    else:
        pdblist.append(pdbid)                
    
# os.chdir("tmp/mol2")
with open("tmp/mol2/listCavTagged", "w") as f:
    for pdbid in pdblist:
        if(pdbid not in blacklisted):
            f.write(f"{pdbid}.mol2\n")
print("computing fingerprints")

ps=subprocess.run(["java","-jar","/home/ubuntu/3decision/pc/fuzcav/FuzCav/dist/3pointPharCav.jar","-d","/home/ubuntu/3decision/pc/fuzcav/FuzCav/utils/resDef/tableDefCA.txt","-t","/home/ubuntu/3decision/pc/fuzcav/FuzCav/utils/triplCav/interval.txt","-l","tmp/mol2/listCavTagged","-o",f"output/{basename}/FPCount.txt","-c"],env={"FuzCav":"/home/ubuntu/3decision/pc/fuzcav/FuzCav"}, check=False, capture_output=True)
print("printing output")

while "NullPointerException" in ps.stderr.decode('utf-8').strip():
    out=ps.stdout.decode('utf-8').strip()
    pdbcode=out[-4:]
    index=pdblist.index(pdbcode)
    if(len(pdblist)>index+1):
        blacklisted.append(pdblist[index+1])
        pdblist.remove(pdblist[index+1])
        # print(pdblist[index+1])

    with open("tmp/mol2/listCavTagged", "w") as f:
        for pdbid in pdblist:
            if(pdbid not in blacklisted):
                f.write(f"{pdbid}.mol2\n")
    print(f"computing fingerprints {len(blacklisted)} vs {len(pdblist)}")

    ps=subprocess.run(["java","-jar","/home/ubuntu/3decision/pc/fuzcav/FuzCav/dist/3pointPharCav.jar","-d","/home/ubuntu/3decision/pc/fuzcav/FuzCav/utils/resDef/tableDefCA.txt","-t","/home/ubuntu/3decision/pc/fuzcav/FuzCav/utils/triplCav/interval.txt","-l","listCavTagged","-o",f"output/{basename}/FPCount.txt","-c"],env={"FuzCav":"/home/ubuntu/3decision/pc/fuzcav/FuzCav"}, check=False, capture_output=True)


with open(f"output/{basename}/complist.txt", "w") as f:
    for pdbid1 in pdblist:
        for pdbid2 in pdblist:
            if(pdbid1 not in blacklisted and pdbid2 not in blacklisted):
                f.write(f"{pdbid1} {pdbid2}\n")

print(os.getcwd())
with open(blacklistname,"w") as f:
    f.writelines([pdbid+"\n" for pdbid in blacklisted])

print("comparing")
with open(outname,"w") as f:
    subprocess.run(["/home/ubuntu/3decision/pc/fuzcav/FuzCav/utils/simC++/simCalc",f"output/{basename}/FPCount.txt",f"output/{basename}/complist.txt"],stdout=f,text=True)
decoydict={}
with open(outname,"r") as f:
    lines=f.readlines()
    for lidx,line in enumerate(lines):
        tmp=line.split("\t")
        # decoy1=str(alldata.loc[alldata["pdbid"]==tmp[0].upper(),"decoy"].values[0])
        # decoy2=str(alldata.loc[alldata["pdbid"]==tmp[1].upper(),"decoy"].values[0])
        pdbid1=tmp[0]
        pdbid2=tmp[1]
        if pdbid1 not in decoydict: 
            decoydict[pdbid1]=str(alldata.loc[alldata["pdbid"]==pdbid1.upper(),"decoy"].values[0])
        if pdbid2 not in decoydict: 
            decoydict[pdbid2]=str(alldata.loc[alldata["pdbid"]==pdbid2.upper(),"decoy"].values[0])
        decoy1=decoydict[pdbid1]
        decoy2=decoydict[pdbid2]
        lines[lidx]=f"{line.strip()}\t{decoy1}\t{decoy2}\n"

with open(outname,"w") as f:
    f.write("pdb_id1\tpdb_id2\tsim\tpdb1decoy\tpdb2decoy\n")
    f.writelines(lines)
import pandas as pd
import numpy as np
import ast
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

# Accession;Residues;SCOP2
# P29274;

    

scop_dict={}
residue_lists=pd.read_csv("residueList.tsv",sep=";")
for index, row in residue_lists.iterrows():
    scop_list=ast.literal_eval(row["SCOP2"])
    if len(scop_list)>=3:
        if scop_list[2] in scop_dict:
            scop_dict[scop_list[2]].append(row["Accession"])
        else:  
            scop_dict[scop_list[2]]=[row["Accession"]]
    else:
        print(f"dropping {row['Accession']} no scop2 classification available")
print(scop_dict)


chembl_data=pd.read_csv("chembl_compounds.csv",sep=";")
chembl_data['mol'] = chembl_data['smiles'].apply(lambda x: Chem.MolFromSmiles(x))
filtered_data = chembl_data.dropna(subset=['mol'])



filtered_data['fp'] = filtered_data['mol'].apply(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 3))

filtered_data = filtered_data.dropna(subset=['fp'])

accessions=filtered_data["accession"].unique()
decoy_candidates={}
final_dict={}
with open("decoy_sets_with_chembl.tsv","w") as f:
    with open("decoy_sets_with_chembl_and_same_SCOP.tsv","w") as f2:
        f.write("Accession;DecoyCandidates\n")
        f2.write("Accession;SimilarAccessions;DecoyCandidates\n")
        for aidx,accession1 in enumerate(accessions):
            candidates=[]
            for accession2 in accessions:
                if(accession1!=accession2):
                    print(accession1,accession2)
                    fps1=list(filtered_data.loc[filtered_data["accession"]==accession1]["fp"])
                    fps2=list(filtered_data.loc[filtered_data["accession"]==accession2]["fp"])
                    pairs=[]
                    for fp1id,fp1 in enumerate(fps1):
                        s = np.where(np.array(DataStructs.BulkTanimotoSimilarity(fp1, fps2)) >= 0.9)
                        if len(s[0]):
                            [pairs.append(el ) for el in s[0]]
                    if(len(pairs)==0):
                        print(pairs)
                        candidates.append(accession2)
                    print(candidates)
            decoy_candidates[accession1]=candidates
            
            decoy_candidate_list = decoy_candidates[accession1]
            tmp=[scop_class for scop_class in scop_dict.keys() if accession1 in scop_dict[scop_class]]
            if len(tmp)>0:                
                scop_class = tmp[0]
                
                print(decoy_candidate_list,scop_class)
                final_decoy_list=[candidate for candidate in decoy_candidate_list if candidate not in scop_dict[scop_class]]
                f.write(f"{accession1};{str(final_decoy_list)}\n")
                f2.write(f"{accession1};{str(scop_dict[scop_class])};{str(final_decoy_list)}\n")
            else:
                print(f"{accession1} without scop class")
            f.flush()
            f2.flush()

        for aidx,accession in enumerate(accessions):
            decoy_candidate_list = decoy_candidates[accession]
            tmp=[scop_class for scop_class in scop_dict.keys() if accession in scop_dict[scop_class]]
            if len(tmp)>0:                
                scop_class = tmp[0]
                
                print(decoy_candidate_list,scop_class)
                final_decoy_list=[candidate for candidate in decoy_candidate_list if candidate not in scop_dict[scop_class]]
                f.write(f"{accession};{str(final_decoy_list)}\n")
                f2.write(f"{accession};{str(scop_dict[scop_class])};{str(final_decoy_list)}\n")
                scop_dict[scop_class].remove(accession)
                final_dict[accession]=(scop_dict[scop_class],final_decoy_list)

with open("decoy_sets_with_chembl_and_same_SCOP_reduced_redundancy.tsv","w") as f3:
    f3.write("Accession;SimilarAccessions;DecoyCandidates\n")
    for accession in final_dict.keys():
        decoys=[]
        forbidden=[]
        for decoyaccession in final_dict[accession][1]:
            if decoyaccession in final_dict.keys():
                forbidden+=final_dict[decoyaccession][0]
            
                if decoyaccession not in forbidden:
                    decoys.append(decoyaccession)
        f3.write(f"{accession};{str(final_dict[accession][0])};{str(decoys)}\n")
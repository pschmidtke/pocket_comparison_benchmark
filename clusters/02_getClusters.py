import json
import requests
from biotite.sequence.io import fasta
import gemmi
from gemmi import cif
# import urllib
import numpy as np
from scipy.spatial.distance import cdist
# import scipy.spatial.distance as ssd
import scipy
# import matplotlib.pyplot as plt
import pandas as pd
# from IPython.display import Markdown
# from tabulate import tabulate
import blosum
subst_matrix = blosum.BLOSUM(62)
from ast import literal_eval
import pickle
import os
import time
import matplotlib.pyplot as plt
import sys

def getRcsbPolymerFromSequence(sequence):
    query="""{
    "query": {
        "type": "terminal",
        "service": "sequence",
        "parameters": {
        "evalue_cutoff": 1e-30,
        "identity_cutoff": 0.1,
        "sequence_type": "protein",
        "value": \""""+sequence+"""\"
        }
    },
    "request_options": {
        "results_verbosity": "verbose",
        "scoring_strategy": "sequence",
        "return_all_hits": true
    },
    "return_type": "polymer_entity"
    }"""

    url="https://search.rcsb.org/rcsbsearch/v2/query?json="+query
    response=requests.get(url)
    hits=json.loads(response.text)["result_set"]
    return(hits)
    # print(hits[0])


def getAuthorBasedChainNames(hits):
    #get all polymer identifiers from previous results
    polymer_identifiers=[hit["identifier"] for hit in hits]

    #now let's build another graphql query against the rcsb to get the author defined chain names of the structures that correspond to the polymer identified. 
    query="""
    {
    polymer_entities(entity_ids: """+str(polymer_identifiers).replace("'","\"")+""")
    {
        rcsb_id
        rcsb_polymer_entity_container_identifiers {
        auth_asym_ids
        }
        rcsb_polymer_entity_align{
        reference_database_name
                aligned_regions{
                entity_beg_seq_id
                length
                ref_beg_seq_id
                }
        }
    }
    }
    """
    # print(query)
    url=f"https://data.rcsb.org/graphql?query={query}"
    try:
      response=requests.get(url)
    except Exception as err:
      print(err)
      time.sleep(120)
      response=requests.get(url)
    chains_tmp=response.json()["data"]["polymer_entities"]

    #here we know what polymer id corresponds to what chain name
    chain_polymer_mapping={chain["rcsb_id"]:chain["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] for chain in chains_tmp}
    # for chain in chains_tmp:
    #     print(chain["rcsb_polymer_entity_align"])
    #     print(len(chain["rcsb_polymer_entity_align"][0]["aligned_regions"]))
    start_residue_mapping={chain["rcsb_id"]:chain["rcsb_polymer_entity_align"][0]["aligned_regions"][0]["ref_beg_seq_id"] for chain in chains_tmp if chain["rcsb_polymer_entity_align"]!=None}
    #From our blast hits, let's keep only the sequence matching bits to simplify the object a bit:
    clean_hits={hit["identifier"]:hit["services"][0]["nodes"][0]["match_context"][0] for hit in hits}
    return(chain_polymer_mapping,start_residue_mapping,clean_hits)

def transformResidueNumbers(residue_numbers, rcsbpolymerhit):
  df = pd.DataFrame(columns = ["aa_query", "aa_subject", "query_residue_number","subject_residue_number"])

  currentindex=0
  subject_aligned_seq_list=list(rcsbpolymerhit["subject_aligned_seq"])
  query_current_residue_number=rcsbpolymerhit["query_beg"]
  subject_current_residue_number=rcsbpolymerhit["subject_beg"]
  
  for idx,queryaa in enumerate(list(rcsbpolymerhit["query_aligned_seq"])):
    if(queryaa!='-'):
      query_residue_number=query_current_residue_number
      query_current_residue_number+=1
    else: 
      query_residue_number=-1
    if(subject_aligned_seq_list[idx]!='-'):
      subject_residue_number=subject_current_residue_number
      subject_current_residue_number+=1
    else:
      subject_residue_number=-1
    new_row=pd.DataFrame({"aa_query":queryaa,"aa_subject":subject_aligned_seq_list[idx],"query_residue_number":query_residue_number,"subject_residue_number":subject_residue_number}, index=[0])
    df=pd.concat([df.loc[:],new_row]).reset_index(drop=True)
  
  subject_residue_list=[df.loc[df['query_residue_number'] == residue_number, 'aa_subject'].values[0] for residue_number in residue_numbers]
  return((df.loc[df['query_residue_number'].isin(residue_numbers),"subject_residue_number"].tolist(),subject_residue_list))


def clusterMatrices(matrixList,selectedResidueList,blosumWeight=1.0):  
  nMatrices=len(matrixList)
  result=np.zeros((nMatrices,nMatrices))
  for i in range(nMatrices):
    for j in range(nMatrices):
      if i==j:
        result[i][j]=0.0
      elif i<j:
        r1=selectedResidueList[i]
        r2=selectedResidueList[j]
        blosum_score=(np.array([subst_matrix[r1[idx]][r2[idx]] for idx in range(len(r1))])<0.0).sum()
        result[i][j]=np.mean(np.abs(matrixList[i]-matrixList[j]))+blosumWeight*blosum_score
        result[j][i]=result[i][j]
  distArray = scipy.spatial.distance.squareform(result)
  
  clusters=scipy.cluster.hierarchy.linkage(distArray, method='single', metric='euclidean')
  corrcoef=scipy.cluster.hierarchy.cophenet(clusters,distArray)
  print("correlation coeff single euclidean")
  sys.stderr.write(f"{corrcoef[0]}\n")
  
  clusters=scipy.cluster.hierarchy.linkage(distArray, method='complete', metric='euclidean')
  corrcoef=scipy.cluster.hierarchy.cophenet(clusters,distArray)
  print("correlation coeff complete euclidean")
  sys.stderr.write(f"{corrcoef[0]}\n")
  
  clusters=scipy.cluster.hierarchy.linkage(distArray, method='ward', metric='euclidean')
  corrcoef=scipy.cluster.hierarchy.cophenet(clusters,distArray)
  print("correlation coeff ward euclidean")
  sys.stderr.write(f"{corrcoef[0]}\n")
  
  clusters=scipy.cluster.hierarchy.linkage(distArray, method='centroid', metric='euclidean')
  corrcoef=scipy.cluster.hierarchy.cophenet(clusters,distArray)
  print("correlation coeff centroid euclidean")
  sys.stderr.write(f"{corrcoef[0]}\n")
  
  clusters=scipy.cluster.hierarchy.linkage(distArray, method='average', metric='euclidean')
  corrcoef=scipy.cluster.hierarchy.cophenet(clusters,distArray)
  print("correlation coeff average euclidean")
  sys.stderr.write(f"{corrcoef[0]}\n")

  return(clusters)



def getContactMatrix(pdbCode, chainCode, residueSelection, debug=False):
  cifpath="/mnt/share/data2/pdb/structures_cifs/"+pdbCode.lower()[1:3]+"/"+pdbCode.lower()+".cif.gz"
  block = cif.read(cifpath)[0]
  structure=gemmi.make_structure_from_block(block)
  positions=[]

  # debug=True
  for model in structure:
    if model.name=="1":
      for chain in model:
        if chain.name == chainCode:
          for residue in chain:
            if residue.label_seq in residueSelection:
              # print(residueSelection)
              if debug:
                print(residue.seqid)
                print(residue)
              for atom in residue:
                if atom.name=="CA":
                  if debug: print("ok")
                  positions.append(atom.pos.tolist())
                  break #we need that for multiple occurences

  if(len(positions)!=len(residueSelection)):
    print("Not all positions found for "+pdbCode+" discarding structure")
    return [None]

  positions_np=np.array(positions)
  return cdist(positions_np, positions_np, 'euclidean')



fasta_file = fasta.FastaFile.read("sequences.fa")
residueLists=pd.read_csv("residueList.tsv",sep=";",header=0)
residueListDict=residueLists.set_index('Accession').T.to_dict('list')

for name,seq in fasta_file.items():
      accession=name.split("|")[1]
      print(accession)
      sys.stderr.write(f"{accession}\n")

      # if accession=="P29274":
      if(not os.path.exists("data/"+accession+"contactMatrices.pkl")):

        # print(accession)
        # if accession=='P07550':
            #get all PDB structures containing this sequence or similar
          hits=getRcsbPolymerFromSequence(seq)
          # hitselected=[hit for hit in hits if hit["identifier"]=="3BBW_1"]
          # hits=hitselected
          (chain_polymer_mapping,start_residue_mapping,clean_hits)=getAuthorBasedChainNames(hits)
          residue_list=np.array(literal_eval(residueListDict[accession][0]),dtype="int")
          # create the contact matrices (CA based for now)
          # print(residue_list)
          contactMatrices=[]
          resultResidueNames=[]
          selectedResidues=[]
          pdbkeys=list(chain_polymer_mapping.keys())
      
          for pdbid in pdbkeys:
            # if pdbid=='6GPX_1':
              appendflag=0
              chain=chain_polymer_mapping[pdbid][0] #take only the first chain
              pdbcode=pdbid.split("_")[0]
              if(residue_list[0]>=clean_hits[pdbid]["query_beg"] and residue_list[-1]<=clean_hits[pdbid]["query_end"]):
                tmp_result=transformResidueNumbers(residue_list,clean_hits[pdbid])
                mapped_residues=np.array(tmp_result[0])
                
                if(os.path.exists("/mnt/share/data2/pdb/structures_cifs/"+pdbcode.lower()[1:3]+"/"+pdbcode.lower()+".cif.gz")):
                  cm=getContactMatrix(pdbcode,chain,mapped_residues)
                  if(len(cm)>1):
                    contactMatrices.append(cm)
                    resultResidueNames.append(tmp_result[1])
                    selectedResidues.append([chain+":"+str(res) for res in mapped_residues])
                    print("mapped residues")
                    print(mapped_residues)
                    appendflag=1
              if appendflag==0:
                contactMatrices.append(None)
                resultResidueNames.append(None)
                selectedResidues.append([None])

          none_indices = [ic for ic, res in enumerate(selectedResidues) if res[0] is None]

          labels=[structure.split("_")[0] for structure in pdbkeys]
          chainCodes=[chain_polymer_mapping[structure][0] for structure in pdbkeys]
          
          

          filteredContactMatrices = [matrix for i, matrix in enumerate(contactMatrices) if i not in none_indices]
          filteredresultResidueNames = [bl for i, bl in enumerate(resultResidueNames) if i not in none_indices]
          filteredSelectedResidues = [bl for i, bl in enumerate(selectedResidues) if i not in none_indices]
          filteredLabels = [label for idx, label in enumerate(labels) if idx not in none_indices]
          filteredChainCodes= [code for idx, code in enumerate(chainCodes) if idx not in none_indices]

          finalList=pd.DataFrame(({"pdbid":filteredLabels,"chain":filteredChainCodes,"residues":filteredSelectedResidues}))
          finalList.to_csv("data/"+accession+"finalList.tsv",sep="\t",index=False)


          file = open('data/'+accession+'contactMatrices.pkl', 'wb')
          pickle.dump({"contactMatrices":filteredContactMatrices,"fileredResidueNames":filteredresultResidueNames,"filteredLabels":filteredLabels,"filteredChainCodes":filteredChainCodes,"finalList":finalList}, file)
          file.close()
          
          clusters=clusterMatrices(filteredContactMatrices,filteredresultResidueNames)
          plt.figure(figsize=(8, len(filteredChainCodes)/4))
          
          scipy.cluster.hierarchy.dendrogram(clusters,labels=filteredLabels,orientation='right',leaf_font_size=15,color_threshold=1.0)
          # ax.tick_params(labelsize=15)
          plt.savefig("data/"+accession+"_dendrogram.png")

      

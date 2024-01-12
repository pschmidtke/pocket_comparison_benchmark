import pandas as pd
import os
import pickle
import numpy as np
import blosum
from scipy.spatial.distance import cdist
subst_matrix = blosum.BLOSUM(62)
import scipy

def getDistanceMatrix(matrixList,selectedResidueList,blosumWeight=1.0):
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
  return(result)


data=pd.read_csv("data.csv",sep=";")
for index, row in data.iterrows():
  # if index==0:
    accession=row["Accession"]
    print(accession)
    clustername=f"data/{accession}contactMatrices.pkl"
    if os.path.exists(clustername):
        with open(clustername,"rb") as f:
            cm=pickle.load(f)
            # print(cm["finalList"])
            # print(cm["contactMatrices"].shape)
            # print(cm)
            dm=getDistanceMatrix(cm["contactMatrices"],cm["fileredResidueNames"])
            clusters=scipy.cluster.hierarchy.linkage(dm, method='single', metric='euclidean')
            # clusters=clusterMatrices(cm["contactMatrices"],cm["fileredResidueNames"])
            flatclusters = scipy.cluster.hierarchy.fcluster(clusters, 2.01, criterion='distance')
            clusterlabels=np.unique(flatclusters)
            allselected=[]
            for cluster_label in clusterlabels:
                clusterpoints=dm[flatclusters==cluster_label]
                
                centroid=np.mean(clusterpoints,axis=0)
                
                distances=cdist(clusterpoints, [centroid])
                selected=cm["finalList"].loc[flatclusters==cluster_label].iloc[np.argmin(distances)]
                
                # print(cm["finalList"].index[cm["finalList"].loc[cm["finalList"]["pdbid"]==selected["pdbid"]]])
                allselected.append(cm["finalList"].index[cm["finalList"]["pdbid"]==selected["pdbid"]].tolist()[0])
            
            print(cm["finalList"].iloc[allselected])
            cm["finalList"].iloc[allselected].to_csv(f"data/{accession}clusterRepresentatives.csv",sep="\t",index=False,header=True)
            
            
            
                

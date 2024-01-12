import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from sklearn.metrics import matthews_corrcoef
from sklearn.metrics import roc_auc_score
from sklearn.metrics import f1_score
from sklearn.metrics import cohen_kappa_score
def custom_sort_key(s):
    parts = s.split('_')
    return float(parts[0]), float(parts[1]), float(parts[2])


if len(sys.argv)!=3:
    sys.exit("USAGE: python analyzeBenchmark.py benchmarkoutput.csv outputimageprefix")
inputfile=sys.argv[1]
outputfile=sys.argv[2]

metrics=[matthews_corrcoef,f1_score,cohen_kappa_score,roc_auc_score]

data=pd.read_csv(inputfile,sep=" ")
data=data.loc[data["pdb_id1"]!=data["pdb_id2"]]
data=data.loc[-((data["pdb1decoy"]==1) & (data["pdb2decoy"]==1))]

for metric in metrics:
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
    thresholds=np.arange(0.0,1.0,0.05)
    scores=[]
    for threshold in thresholds:
        y_pred=data["sim"]>=threshold
        y_true=(data["pdb1decoy"]==0) & (data["pdb2decoy"]==0)
        if(metric.__name__=="f1_score"):
            score=metric(y_true,y_pred,pos_label=True)
        else:
            score=metric(y_true,y_pred)
        scores.append(score)
    axes.plot(thresholds,scores)


    plt.savefig(f"{outputfile}_{metric.__name__}.png")

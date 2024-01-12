import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from sklearn.metrics import matthews_corrcoef
def custom_sort_key(s):
    parts = s.split('_')
    return float(parts[0]), float(parts[1]), float(parts[2])


if len(sys.argv)!=3:
    sys.exit("USAGE: python analyzeBenchmark.py benchmarkoutput.csv outputimage.png")
inputfile=sys.argv[1]
outputfile=sys.argv[2]

data=pd.read_csv(inputfile,sep="\t")
data=data.loc[data["pdb_id1"]!=data["pdb_id2"]]
data=data.loc[-((data["pdb1decoy"]==1) & (data["pdb2decoy"]==1))]

fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))

thresholds=np.arange(0.0,1.0,0.05)

mccs=[]
for threshold in thresholds:
    y_pred=data["sim"]>=threshold
    y_true=(data["pdb1decoy"]==0) & (data["pdb2decoy"]==0)
    mcc_score=matthews_corrcoef(y_true,y_pred)
    mccs.append(mcc_score)
axes.plot(thresholds,mccs)


plt.savefig(outputfile)

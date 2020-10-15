import numpy as np
import pandas as pd
import pandas 
import json
import argparse
import os

parser=argparse.ArgumentParser()
parser.add_argument("--json",type=str,default=None)
parser.add_argument("--label",type=str,default=None)
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--meta",type=str,default=None)
args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)


colors=["#0080ff","#0028ff","#0000de","#00007f","#36ffc0","#00d4ff","#c0ff36","#7cff79","#ff9400","#ffe500","#de0000","#ff4600"]

print("### Loading json file")
fread=open(args.json,"r")
data=json.load(fread)
df=pd.read_csv(args.label)
label_list=df["label"].values.tolist()

############################
n_label=len(np.unique(label_list))
color=np.random.choice(colors,n_label,replace=False)
label_colors={label:color for label,color in zip(np.unique(label_list),color)}

label_npy=np.array(label_list)
with open(os.path.join(args.outdir,"label.npy"),"wb") as f:
    np.save(f,label_npy)
f.close()

meta=pd.read_csv(args.meta)
sample_list=meta["Sample"].tolist()
status_list=meta["Status"].tolist()

#############################
n=len(np.unique(sample_list))
color=np.random.choice(colors,n,replace=False)
sample_colors={label:color for label,color in zip(np.unique(sample_list),color)}

sample_npy=np.array(sample_list)
with open(os.path.join(args.outdir,"sample.npy"),"wb") as f:
    np.save(f,sample_npy)
f.close()

#############################
n=len(np.unique(status_list))
color=np.random.choice(colors,n,replace=False)
status_colors={label:color for label,color in zip(np.unique(status_list),color)}

status_npy=np.array(status_list)
with open(os.path.join(args.outdir,"status.npy"),"wb") as f:
    np.save(f,status_npy)
f.close()



data["Label"]={"label_colors":label_colors,"label_list":label_list}
data["Sample"]={"label_colors":sample_colors,"label_list":sample_list}
data["Status"]={"label_colors":status_colors,"label_list":status_list}
with open(os.path.join(args.outdir,os.path.basename(args.json)),"w") as f:
    json.dump(data,f)


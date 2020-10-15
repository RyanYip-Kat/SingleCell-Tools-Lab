import scanpy as sc
import numpy as np
import pandas as pd

import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys


from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--files",nargs="+",type=str,default=None)
parser.add_argument("--names",nargs="+",type=str,default=None)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("*** Loading data")
adata=sc.read_h5ad(args.data)
cells=adata.obs_names

files=args.files
names=args.names
assert len(files)==len(names)

for name,f in zip(names,files):
    print("Loadin {} embeding files from {}".format(name,f))
    data=pd.read_csv(f,sep=",")

    assert data.shape[0]==len(cells)
    data.index=data.Barcode.to_list()
    data.drop(["Barcode"],axis=1,inplace=True)
    data=data.loc[cells,:]
    data=np.array(data)
    print(data)

    emb="X_"+name
    print("Adding embedding : {} into adata".format(emb))
    adata.obsm[emb]=data

print("Update Embedding done!")
print("Save...")
adata.write(os.path.join(args.outdir,"adata.h5ad"))



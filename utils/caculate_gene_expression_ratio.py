import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None)
parser.add_argument("-g","--gene",type=str,default="BSG")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)


gene=args.gene
adata=sc.read_h5ad(args.data)
DATA=adata.copy()

DATA=DATA.raw.to_adata()
df=DATA.to_df()
df=df[[gene]]
df[["idents"]]=DATA.obs[["idents"]]
df[["status"]]=DATA.obs[["status"]]

status=np.unique(df[["status"]])
idents=np.unique(df[["idents"]])
for s in status:
    X=df[df["status"].isin([s])]
    N=X.shape[0]
    writer=open(os.path.join(args.outdir,s+"_table.txt"),"w")
    for i in idents:
        y=X[X["idents"].isin([str(i)])]
        m=y.shape[0]
        n=y[y[gene]>0].shape[0]
        
        line="ident : {} expressed : {}/{} in {} [ {} ]".format(i,n,m,s,N)
        writer.write(line+'\n')
    writer.close()


import scanpy as sc
import pandas as pd
import numpy as np
import os
import argparse
import scipy
from scipy import stats

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

if len(status)==1:
    raise ValueError("status length must be > = 2")
elif len(status)==2:
    x=df[df["status"].isin([status[0]])][gene].values
    y=df[df["status"].isin([status[1]])][gene].values
    diff,pvals=stats.ttest_ind(x,y)
    print("The mean difference between : {} and {} is : {},and the pvals is {}".format(status[0],status[1],diff,pvals))

def gene_comparison(adata,gene,by="status"):
    assert by in ["status","idents"]
    DATA=adata.copy()
    DATA=DATA.raw.to_adata()
    df=DATA.to_df()
    df=df[[gene]]
    df[[by]]=DATA.obs[[by]]
    status=np.unique(df[[by]])

    if len(status)==1:
        raise ValueError("status length must be > = 2")
    elif len(status)==2:
        func=stats.ttest_ind
        x=df[df["status"].isin([status[0]])][gene].values
        y=df[df["status"].isin([status[1]])][gene].values
        diff,pvals=func(x,y)
        print("t-test,The mean difference between : {} and {} is : {},and the pvals is {}".format(status[0],status[1],diff,pvals))
    elif len(status)>2:
        func=stats.f_oneway
        v_list=[]
        for s in status:
            x=df[df["status"].isin([s])][gene].values
            v_list.append(x)

        diff,pvals=func(v_list)
        print("oneway anova,The mean difference across status is :{},and the pvals is {}".format(diff,pvals))

    return diff,pvals

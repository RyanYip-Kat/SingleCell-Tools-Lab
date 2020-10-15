import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)

parser.add_argument("--n_top",type=int,default=None,help="the top genes return according to FC")
args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("# Loading data")
adata=sc.read_h5ad(args.data)

def get_markers(adata,n_top=None,pval=None,logfc=None):
    varms=adata.varm["de_res"]
    markers_df=pd.DataFrame(varms)
    columns=markers_df.columns

    fields=np.unique([column.split(":")[0] for column in columns])
    groups=np.unique([column.split(":")[1] for column in columns])
    names=adata.var_names.tolist()
    table=[]
    for group in groups:
        try:
           print("Get markers from : {}".format(group))
           col=[field+":"+str(group) for field in fields]
           marker_df=markers_df[col]
           marker_df.columns=fields
           marker_df["cluster"]=group
           marker_df["names"]=names
           marker_df=marker_df.sort_values("log_fold_change",ascending=False)
           if pval is not None:
               marker_df=marker_df[marker_df.t_qval < pval]
           if logfc is not None:
               marker_df=marker_df[marker_df.log_fold_change > logfc]
           if n_top is not None:
              marker_df=marker_df[:n_top]

           table.append(marker_df)
        except:
           print("group : {} has no DE".format(group))
           
    mat=pd.concat(table,axis=0)
    if n_top is not None:
       filename=os.path.join(args.outdir,str(n_top)+"_markers.csv")
    else:
       filename=os.path.join(args.outdir,"_markers.csv")
    
    print("Save markers in {}".format(filename))
    mat.to_csv(filename,sep=",",index=False)



get_markers(adata,n_top=args.n_top,pval=0.05,logfc=0.25)
get_markers(adata,n_top=500,pval=0.05,logfc=0.25)

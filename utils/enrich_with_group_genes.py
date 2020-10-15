import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default="leiden",help="data")
parser.add_argument("--column2",type=str,default=None,help="data")
parser.add_argument("--groupby",type=str,default="status",help="data")
parser.add_argument("--method",type=str,default="t-test",help="t-test_overestim_var,wilcoxon")
parser.add_argument("--n_genes",type=int,default="50")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

org="hsapiens"
adata=sc.read_h5ad(args.data)
if args.subset1 is not None and args.column1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)

    #adata=adata.raw.to_adata()
    #adata.raw=adata
    sc.pp.normalize_total(adata, target_sum=1e4)

sc.tl.rank_genes_groups(adata,groupby=args.groupby,method=args.method,n_genes=args.n_genes,corr_method="bonferroni")
groups=np.unique(adata.obs[args.groupby])
tables=[]

print("Enrich")
for g in groups:
    t=sc.get.rank_genes_groups_df(adata,group=g)
    genes=t.names.tolist()
    df=sc.queries.enrich(genes,org=org)
    df["cluster"]=g
    tables.append(df)

m=pd.concat(tables,axis=0)
filename=os.path.join(args.outdir,args.groupby+"_enrich.csv")
m.to_csv(filename,sep=",",index=False)

print("Done")



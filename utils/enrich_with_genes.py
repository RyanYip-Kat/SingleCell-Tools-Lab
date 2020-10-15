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
parser.add_argument("--org",type=str,default="hsapiens")
#parser.add_argument("--genes",nargs="+",type=str,default=None)

parser.add_argument("--genes",type=str,default=None)
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

df=pd.read_csv(args.genes,sep="\t",header=None)
marker_genes=df.loc[:,0].values.tolist()

df=sc.queries.enrich(marker_genes,org=args.org)
filename=os.path.join(args.outdir,"enrich_genes.csv")
df.to_csv(filename,sep=",",index=False)

print("Done")



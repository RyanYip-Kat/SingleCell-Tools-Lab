import argparse
import os
import scanpy as sc
import pandas as pd


parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--column",type=str,default="leiden")



args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("*** Loading data")
adata=sc.read_h5ad(args.data)
metadata=adata.obs

if args.column in metadata.columns:
    if "Sample" in metadata.columns:
        df=pd.crosstab(metadata["Sample"],metadata[args.column])
        df.to_csv(os.path.join(args.outdir,"sample_"+args.column+"_crosstable.csv"),sep=",")

    if "Status" in metadata.columns:
        df=pd.crosstab(metadata["Status"],metadata[args.column])
        df.to_csv(os.path.join(args.outdir,"status_"+args.column+"_crosstable.csv"),sep=",")

else:
    raise ValueError("Please input valid column in metadata")

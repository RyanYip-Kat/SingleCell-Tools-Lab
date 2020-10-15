import scanpy as sc
import os
import pandas as pd
import argparse
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--barcode",type=str,default=None,help="the csv file of barcode to be specified")
parser.add_argument("--column",type=str,default="celltype",help="the columns in adata.obs to be used")
parser.add_argument("--name",type=str,default=None)
parser.add_argument("--key",type=str,default="clusters")

args=parser.parse_args()

if args.outdir is not None:
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

print("*** Loading data from : {}".format(args.data))
adata=sc.read_h5ad(args.data)
DATA=pd.read_csv(args.barcode,sep=",")

barcodes=DATA["Barcode"].values.tolist()
metadata=adata.obs
metadata["cells"]=metadata.index

print(" ")
print("---------")
print("### Before add counter")
print(Counter(metadata[args.column].values.tolist()))
new_clusters=[]
for barcode,celltype in zip(metadata["cells"].values.tolist(),metadata[args.column].values.tolist()):
    if barcode in barcodes:
        new_clusters.append(args.name)
    else:
        new_clusters.append(celltype)

print(" ")
print("---------")
print("### After add counter")
print(Counter(new_clusters))
adata.obs[args.key]=new_clusters
if args.outdir is not None:
    adata.write(os.path.join(args.outdir,"adata.h5ad"))
else:
    adata.write(args.data)



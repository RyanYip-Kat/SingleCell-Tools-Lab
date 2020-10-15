import argparse
import os
import scanpy as sc
import numpy as np
#from core.colors import *  # color dict

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default=None)
parser.add_argument("-c","--column",type=str,default="leiden")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("*** Loading data")
adata=sc.read_h5ad(args.data)
metadata=adata.obs

assert args.column in metadata.columns
cell_label=metadata[args.column]
cell_label.to_csv(os.path.join(args.outdir,"cell_label.tsv"),sep="\t",index=False,header=None)


label_color={"0": "#FF0000",
             "1": "#836FFF",
             "2": "#0000FF",
             "3": "#C6E2FF",
             "4": "#548B54",
             "5": "#00FF00",
             "6": "#FFF68F",
             "7": "#8B864E",
             "8": "#FFFF00",
             "9": "#FFD700",
             "10": "#8B658B",
             "11": "#FF6A6A",
             "12": "#FFD39B",
             "13": "#EE2C2C",
             "14": "#BF3EFF",
             "15": "#483D8B",
             "16": "#0000CD",
             "17": "#00CED1",
             "18": "#9ACD32",
             "19": "#EEB422",
             "20": "#8B8989",
             "21": "#008B8B"
             }

adata.obs["label_color"] = [label_color[x] for x in adata.obs[args.column]]
cell_label_color=adata.obs[[args.column,"label_color"]]
cell_label_color.to_csv(os.path.join(args.outdir,"cell_label_color.tsv"),sep="\t",index=False,header=None)

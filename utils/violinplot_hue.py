import scanpy as sc
import argparse
import sys
import os
import seaborn as sns
sys.path.append("../")

from core.plotting import violin_hue

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-g","--genes",nargs="+",type=str,default=None)
parser.add_argument("-G","--groupby",type=str,default="celltype")
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

genes=args.genes
adata=sc.read_h5ad(args.data)

sns.set(style="whitegrid", palette="pastel", color_codes=True)
for gene in genes:
    violin_hue(adata,groupby=args.groupby,
            gene_name=gene,
            hue="status",
            hue_order=["YA","AA"],
            use_raw=True,
            figdir=args.outdir)

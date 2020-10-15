import scanpy as sc
import argparse
import sys
import os
sys.path.append("../")

from core.utils import get_rank_group_genes

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default="output")
parser.add_argument("-d","--data",type=str,default=None,help="data")
parser.add_argument("-n","--n_top",type=int,default=50)
parser.add_argument("-k","--key",type=str,default="rank_genes_groups")
args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
get_rank_group_genes(adata,outdir=args.outdir,n_top=args.n_top,key=args.key)
get_rank_group_genes(adata,outdir=args.outdir,key=args.key)

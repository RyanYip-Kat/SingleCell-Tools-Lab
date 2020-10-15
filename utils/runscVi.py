import argparse
import os
import scanpy as sc
import scanpy.external as sce
import numpy as np
import sys

sys.path.append("../")

from core.utils import subset_by_column,subset_by_cell_feature

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--batch_key",type=str,default=None)
parser.add_argument("--n_epochs",type=int,default=400)
parser.add_argument("--lr",type=float,default=0.001)

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
sc.settings.n_jobs=8
print("*** Loading data")
adata=sc.read_h5ad(args.data)
adata=adata.raw.to_adata()
if args.batch_key is not None:
    assert args.batch_key in adata.obs.columns

print("Training with scvi...")
sce.pp.scvi(adata,n_hidden=128,
        n_latent=10,
        n_layers=1,
        n_epochs=args.n_epochs,
        lr=args.lr,
        batch_key=args.batch_key,
        use_cuda=False)

print("Save..")
adata.write(os.path.join(args.outdir,"adata.h5ad"))



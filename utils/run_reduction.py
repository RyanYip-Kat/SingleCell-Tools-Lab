import argparse
import os
import scanpy as sc
import numpy as np
import sys

sys.path.append("../")
from core.functions import batch_correct
parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--use_rep",type=str,default="X_harmony",choices=["X_pca","X_harmony","X_scvi"],help="The embedding for tsne,umap,clusters")
parser.add_argument("--resolution",type=float,default=1.2)
parser.add_argument("--batch", action='store_true', default=False,help="batch correct")
parser.add_argument("--batch_method",type=str,default="combat",choices=["combat","mnn","bbknn"])
parser.add_argument("--batch_key",type=str,default="idents")
args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
sc.settings.n_jobs=12
print("# Loading data")
adata=sc.read_h5ad(args.data)

def reduction(adata,resolution=1.5,use_rep="X_hamony",batch_method="combat",batch_key="idents",batch=True):

      # use_rep : X_pca,X_harmony,X_scvi...
      assert use_rep in adata.obsm_keys()
      if batch:
          adata=batch_correct(adata,batch_key=batch_key,method=batch_method)

      sc.tl.pca(adata,n_comps=50,svd_solver='arpack')
      print("# Run neighborhood graph ")
      sc.pp.neighbors(adata, n_neighbors=20,use_rep=use_rep)

      print("# Run tsne")
      sc.tl.tsne(adata,use_rep=use_rep,learning_rate=200)

      print("# Clustering")
      sc.tl.leiden(adata,resolution=resolution)
      #sc.tl.louvain(adata,resolution=resolution)

      print("# Run umap")
      sc.tl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
      sc.pl.paga(adata, plot=False)
      sc.tl.umap(adata, init_pos='paga')

      return adata

adata=reduction(adata,resolution=args.resolution,use_rep=args.use_rep,batch_method=args.batch_method,batch_key=args.batch_key,batch=args.batch)
sc.tl.rank_genes_groups(adata,groupby="leiden", method='t-test',n_genes=200,corr_method="bonferroni")
filename=os.path.join(args.outdir,"adata.h5ad")
print("# Save adata into : {}".format(filename))
adata.write(filename)

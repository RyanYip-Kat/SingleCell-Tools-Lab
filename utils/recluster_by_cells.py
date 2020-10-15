import argparse
import os
import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import sys

sys.path.append("../")

from core.tl import reduction,process
from core.core import Model
from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--path",type=str,default=None)
parser.add_argument("--cells",type=str,default=None)
parser.add_argument("--batch_key",type=str,default=None)
parser.add_argument("--invert", action='store_true', default=False)
parser.add_argument("--batch_correct", action='store_true', default=False,help="batch correct")
parser.add_argument("--batch_method",type=str,default="combat",choices=["combat","mnn","bbknn","harmony"])
parser.add_argument("--resolution",type=float,default=1.2)
parser.add_argument("--use_rep",type=str,default="X_pca")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

########################
def batch_correct(adata,batch_key="idents",method="combat"):
    #print("# correct batch effect")
    #sc.pp.normalize_total(adata, target_sum=1e4)
    #sc.pp.log1p(adata)
    #sc.pp.highly_variable_genes(adata,n_top_genes=2000)
    #sc.pp.scale(adata, max_value=10)

    #adata.raw = adata
    if "X_pca" not in adata.obsm_keys():
        sc.tl.pca(adata,n_comps=50,svd_solver='arpack')
    if batch_key is not None:
        assert batch_key in adata.obs.columns
        if method=="combat":
            sc.pp.combat(adata,key=batch_key)
        elif method=="mnn":
            highly_variable_genes = adata.var["highly_variable"]
            hvg=highly_variable_genes.index.to_list()
            adata=sce.pp.mnn_correct(adata,batch_key=batch_key,var_subset=hvg)[0][0]
        elif  method=="bbknn":
            sce.pp.bbknn(adata,batch_key=batch_key,approx=True)
        elif method=="harmony":
            sce.pp.harmony_integrate(adata,key=batch_key,basis="X_pca",adjusted_basis="X_harmony")
        else:
            raise ValueError("Invalid batch effect correct method input!!!")
    else:
        print("# Don't correct batch effect")
    return adata

sc.settings.figdir=args.outdir
sc.settings.n_jobs=8
print("*** Loading data")
if args.path is not None:
    adata=load_data(args.path)
    adata,adata_original=preprocess(adata)
else:
    adata=sc.read_h5ad(args.data)

print("*** subset")
df=pd.read_csv(args.cells)
cells=df.iloc[:,0].values.tolist()
adata=subset_by_cell_feature(adata,cells=cells,invert=args.invert)
adata=adata.raw.to_adata()

if args.batch_correct and args.batch_key is not None:
    adata=batch_correct(adata,args.batch_key,args.batch_method)

print("*** tranform and reduction")
#adata=process(adata,args.batch_key)
adata=reduction(adata,resolution=args.resolution,use_rep=args.use_rep)
print(Counter(adata.obs["leiden"]))

print("*** find markers")
sc.tl.rank_genes_groups(adata,groupby="leiden", method='t-test',n_genes=500,corr_method="bonferroni")
adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')

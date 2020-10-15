import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys
import scanpy.external as sce

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--path",type=str,default=None)
parser.add_argument("--meta",type=str,default=None,help="TGS_union_matrix/min_cells_3/raw/metadata.csv")
parser.add_argument("--batch_key",type=str,default=None)
parser.add_argument("--method",type=str,default="combat",choices=["combat","mnn","bbknn","harmony"])
parser.add_argument("--subset1",nargs="+",type=str,default=None)
parser.add_argument("--column1",type=str,default="leiden")
parser.add_argument("--subset2",nargs="+",type=str,default=None)
parser.add_argument("--column2",type=str,default="leiden")
parser.add_argument("--subset3",nargs="+",type=str,default=None)
parser.add_argument("--column3",type=str,default="leiden")
parser.add_argument("--n_top",type=int,default=2000)
parser.add_argument("--min_cells",type=int,default=200)
args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)
########################
def subset_by_column(adata,subset,col_use="orig_cluster",invert=False):
      metadata=adata.obs
      assert col_use in metadata.columns
      adata.obs[col_use].astype("category")
      subset=[str(s) for s in subset]
      cols=np.unique(adata.obs[col_use].values.tolist())
      print(subset)
      for s in subset:
          if s not in cols:
              raise ValueError("Invalid value in column selected")

      if invert:
          adata_subset=adata[~adata.obs[col_use].isin(subset)]
      else:
          adata_subset=adata[adata.obs[col_use].isin(subset)]
      print("after subset,adata shape is [{},{}]".format(adata_subset.shape[0],adata_subset.shape[1]))
      return adata_subset

def batch_correct(adata,batch_key="idents",method="combat",n_top=2000,use_hvg=False):
    print("# correct batch effect")
    sc.pp.highly_variable_genes(adata,n_top_genes=n_top)
    if use_hvg:
        adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata, max_value=10)

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

########################
print("### Loading dataset")
adata=sc.read_10x_mtx(args.path)
print(adata.shape)
sc.pp.filter_genes(adata,min_cells=args.min_cells)
print(adata.shape)
meta=pd.read_csv(args.meta)
meta.index=meta["cell_id"].values.tolist()
meta=meta.loc[adata.obs_names,:]
adata.obs=meta

time_point_dict={"B0h1":"0h","B0h2":"0h","B3h1":"3h","B3h2":"3h","B6h1":"6h","B6h2":"6h"}
time_point=[time_point_dict[c] for c in adata.obs["Sample"]]
adata.obs["time_point"]=time_point

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

if args.column1 is not None and args.subset1 is not None:
    print("### Subset data by {}".format(args.column1))
    adata=subset_by_column(adata,args.subset1,args.column1,invert=False)
    if args.column2 is not None and args.subset2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2,invert=False)
        if args.column3 is not None and args.subset3 is not None:
            adata=subset_by_column(adata,args.subset3,args.column3,invert=False)

########################

if args.method and args.batch_key is not None:
    adata=batch_correct(adata,args.batch_key,args.method,args.n_top,True)

matrix=adata.to_df().transpose()
print("### Export adata to matrix dataframe")
matrix.to_csv(os.path.join(args.outdir,"matrix.csv"),sep=",")

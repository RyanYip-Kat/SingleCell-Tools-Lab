import scanpy as sc
import pandas as pd
import numpy as np

import argparse
import os
import sys

############################ functions

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

def subset_by_cell_feature(adata,cells=None,features=None,invert=False):
    data=adata.copy()
    data.obs["cells"]=data.obs_names
    Features=adata.var_names
    print("before subset,there are : {} #cells".format(data.shape[0]))
    if cells is not None:
        if invert:
            data=data[~data.obs["cells"].isin(cells)]
        else:
            data=data[data.obs["cells"].isin(cells)]
    if features is not None:
       if invert:
           feature_idx=[False if feature in features else True for feature in Features]
       else:
           feature_idx=[True if feature in features else False for feature in Features]
       data=data[:,feature_idx]
    print("after subset,there are : {} #cells".format(data.shape[0]))
    return data


######################################


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=str,default="output")
    parser.add_argument("--data",type=str,default=None,help="cellranger count and aggr output ")

    parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset clusters in column for downstream")
    parser.add_argument("--column",type=str,default=None,help="the column in metadata using for subset,eg:orig_cluster")
    parser.add_argument("--reverse", action='store_true', default=False,help="reverse for subset if true")

    parser.add_argument("--groupby",type=str,default="leiden_labels",help="data")
    parser.add_argument("--method",type=str,default="t-test",help="t-test_overestim_var,wilcoxon")
    parser.add_argument("--n_genes",type=int,default=500,help="the number genes return")
    args=parser.parse_args()

    if not os.path.exists(args.outdir):
       os.makedirs(args.outdir)
    
    sc.settings.autoshow=False
    sc.settings.figdir=args.outdir
    
    print("# Loading data")
    adata=sc.read_h5ad(args.data)
    if args.subset is not None and args.column is not None:
        adata=subset_by_column(adata,args.subset,args.column,args.reverse)
        sc.pp.normalize_total(adata, target_sum=1e4)
    print("# Run rank genes with group : {},method :{},and return number :{} genes".format(args.groupby,args.method,args.n_genes))
    sc.tl.rank_genes_groups(adata,groupby=args.groupby,method=args.method,n_genes=args.n_genes,corr_method="bonferroni")
    print("# save")
    adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')

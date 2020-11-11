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


def load_data(path,column=None,subset=None,reverse=False):
      graphclust_path=os.path.join(path,"outs/analysis/clustering/graphclust/clusters.csv")
      km_path=os.path.join(path,"outs/analysis/clustering/kmeans_10_clusters/clusters.csv")
      mtx_path=os.path.join(path,"outs/filtered_feature_bc_matrix")

      adata=sc.read_10x_mtx(mtx_path,cache=True)
      orig_cluster=pd.read_csv(graphclust_path)
      km_cluster=pd.read_csv(km_path)

      print("# Add original message")
      adata.obs["orig_cluster"]=[str(i) for i in orig_cluster.Cluster.to_list()]
      adata.obs["km_cluster"]=[str(i) for i in km_cluster.Cluster.to_list()]

      cells=adata.obs_names.to_list()
      idents=[cell.split("-")[1] for cell in cells]
      adata.obs["idents"]=idents

      n_cells=adata.shape[0]
      n_features=adata.shape[1]
      print("# Origin adata shape is [{},{}]".format(n_cells,n_features))

      if column is not None and subset is not None:
          assert column in adata.obs.columns
          adata=subset_by_column(adata,subset,column,reverse)
          print("# After subset adata shape is [{},{}]".format(adata.shape[0],adata.shape[1]))

      adata.obs['n_counts'] =adata.X.sum(axis=1).A1
      adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1
      adata.var['mt'] = adata.var_names.str.startswith('MT-')
      adata.var['rpl'] = adata.var_names.str.startswith('RPL')
      adata.var['rps'] = adata.var_names.str.startswith('RPS')
      sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)
      sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=True,multi_panel=True,size=0.25,show=False,save="_QC.pdf",figsize=[16,10])

      print("# Remove rpl,rps,point genes")
      rpl_genes = adata.var_names.str.startswith('RPL')
      adata=adata[:,~rpl_genes]

      rps_genes = adata.var_names.str.startswith('RPS')
      adata=adata[:,~rps_genes]

      point_genes=np.array([True if "." in var else False for var in adata.var_names])
      adata=adata[:,~point_genes]

      aggregation=os.path.join(path,"outs/aggregation.csv")
      aggr=pd.read_csv(aggregation,sep=",")
      aggr["idents"]=list(range(1,aggr.shape[0]+1))
      status_dict={str(i):s for i,s in zip(aggr["idents"].values.tolist(),aggr["status"].values.tolist())}

      status=[status_dict[i] for i in  adata.obs["idents"]]
      adata.obs["status"]=status
      return adata
######################################


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=str,default="output")
    parser.add_argument("--path",type=str,default=None,help="cellranger count and aggr output ")

    parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset clusters in column for downstream")
    parser.add_argument("--column",type=str,default=None,help="the column in metadata using for subset,eg:orig_cluster")

    parser.add_argument("--reverse", action='store_true', default=False,help="reverse for subset if true")
    args=parser.parse_args()

    if not os.path.exists(args.outdir):
       os.makedirs(args.outdir)
    
    sc.settings.autoshow=False
    sc.settings.figdir=args.outdir
    adata=load_data(args.path,args.column,args.subset,args.reverse)
    print("### save")
    adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')

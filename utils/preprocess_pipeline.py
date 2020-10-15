import scanpy as sc
import pandas as pd
import numpy as np

import argparse
import os

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--path",type=str,default=None,help="cellranger count and aggr output ")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.autoshow=False
sc.settings.figdir=args.outdir


mtx_path=os.path.join(args.path,"outs/filtered_feature_bc_matrix")
graphclust_path=os.path.join(args.path,"outs/analysis/clustering/graphclust/clusters.csv")
km_path=os.path.join(args.path,"outs/analysis/clustering/kmeans_10_clusters/clusters.csv")

print("### Loading data")
adata=adata=sc.read_10x_mtx(mtx_path,cache=True)

orig_cluster=pd.read_csv(graphclust_path)
km_cluster=pd.read_csv(km_path)

print("### Add original message")
adata.obs["orig_cluster"]=[str(i) for i in orig_cluster.Cluster.to_list()]
adata.obs["km_cluster"]=[str(i) for i in km_cluster.Cluster.to_list()]
cells=adata.obs_names.to_list()
idents=[cell.split("-")[1] for cell in cells]
adata.obs["idents"]=idents

print("### Preprocess")
adata.obs['n_counts'] =adata.X.sum(axis=1).A1
adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('mt-')
adata.var['rpl'] = adata.var_names.str.startswith('Rpl')
adata.var['rps'] = adata.var_names.str.startswith('Rps')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)

sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)
      #X=X[(X.obs.n_genes_by_counts> 200) &  (X.obs.n_genes_by_counts <3500),:]
print("Save QC Metric plot")
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=True,multi_panel=True,size=0.25,show=False,save="_QC1.pdf",figsize=[16,10])

sc.pl.violin(adata, ['pct_counts_rpl','pct_counts_rps'],
            jitter=True,multi_panel=True,size=0.25,show=False,save="_QC2.pdf",figsize=[16,10])
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',save="_QC3.pdf",show=False)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',save="_QC4.pdf",show=False)
#print("### Remove dropout features")
#mito_genes = adata.var_names.str.startswith('MT-')
#adata=adata[:,~mito_genes]
#rpl_genes = adata.var_names.str.startswith('RPL')
#adata=adata[:,~rpl_genes]

#rps_genes = adata.var_names.str.startswith('RPS')
#adata=adata[:,~rps_genes]

#print("### data transform")
#sc.pp.normalize_total(adata, target_sum=1e4)
#sc.pp.log1p(adata)
#adata.raw = adata

#sc.pp.highly_variable_genes(adata,n_top_genes=2000)
#sc.pp.regress_out(adata, ['n_counts'])
#adata = adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['n_counts'])

#print("### save")
adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')

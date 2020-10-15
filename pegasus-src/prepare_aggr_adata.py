import scanpy as sc
import pandas as pd
import numpy as np

import argparse
import os
import sys
import re

############################ functions
import scrublet as scr
import scipy.io
import pickle
import matplotlib.pyplot as plt

class DoubletModule(object):
    def __init__(self,
            adata,
            min_counts=2,
            min_cells=3,
            expected_doublet_rate=0.05,
            n_prin_comps=30,
            outdir="./doublet"):

        self._min_counts=min_counts
        self._min_cells=min_cells
        self._expected_doublet_rate=expected_doublet_rate
        self._n_prin_comps=n_prin_comps
        self._outdir=outdir

        self.adata=adata
        if not os.path.exists(self._outdir):
            os.makedirs(self._outdir)

    def detect(self):
        #try:
        #    self.adata=sc.read_h5ad(self._adata)
        #except:
        #    self.adata=sc.read_10x_mtx(self._adata)

        print("### Initialize Scruble")
        counts_matrix=self.adata.raw.X
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=self._expected_doublet_rate)#, sim_doublet_ratio=1.0)

        print("### Detect ,Nomalize,PCA")
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=self._min_counts,
                                                                  min_cells=self._min_cells,
                                                                  min_gene_variability_pctl=85,
                                                                  n_prin_comps=self._n_prin_comps)
        self.adata.obs["doublet_scores_obs"]=doublet_scores
        self.adata.obs["predicted_doublets"]=predicted_doublets
        self.adata.obs["doublet_errors_obs"]=scrub.doublet_errors_obs_
        #self.adata.obs["doublet_errors_sim"]=scrub.doublet_errors_sim_
        #self.adata.obs["doublet_scores_sim"]=scrub.doublet_scores_sim_

        self.scrub=scrub

    def plot_histogram(self,threshold=0.25):
        self.scrub.call_doublets(threshold)
        self.scrub.plot_histogram()
        plt.savefig(os.path.join(self._outdir,"plot_histogram.pdf"),figsize=(12, 12))
        plt.close()

    def set_embedding(self):
        print("#### Run UMAP")
        self.scrub.set_embedding('UMAP', scr.get_umap(self.scrub.manifold_obs_, 30, min_dist=0.3))

        # # Uncomment to run tSNE - slow
        print('#### Running tSNE...')
        self.scrub.set_embedding('tSNE', scr.get_tsne(self.scrub.manifold_obs_, angle=0.9,verbose=True))
        #print('Done.')
        #with open(os.path.join(self._outdir,"scrub.pkl"),"wb") as f:
        #    pickle.dump(self.scrub, f)
        #f.close()

    def add_embedding(self,export=False):
        umap_mat=self.scrub._embeddings["UMAP"]
        tsne_mat=self.scrub._embeddings["tSNE"]

        self.adata.obsm["X_umap_scr"]=umap_mat
        self.adata.obsm["X_tsne_scr"]=tsne_mat
        if export:
            umap_df=pd.DataFrame(umap_mat,index=self.adata.obs_names,columns=["UMAP_1","UMAP_2"])
            tsne_df=pd.DataFrame(tsne_mat,index=self.adata.obs_names,columns=["tSNE_1","tSNE_2"])
            umap_df.to_csv(os.path.join(self._outdir,"umap.csv"),sep=",")
            tsne_df.to_csv(os.path.join(self._outdir,"tsne.csv"),sep=",")
    def save(self):
        self.adata.write(os.path.join(self._outdir,"adata_doublet.h5ad"), compression='gzip')



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
      adata.var['mt'] = adata.var_names.str.startswith('mt-')
      adata.var['rpl'] = adata.var_names.str.startswith('Rpl')
      adata.var['rps'] = adata.var_names.str.startswith('Rps')
      sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)
      sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=True,multi_panel=True,size=0.25,show=False,save="_QC.pdf",figsize=[16,10])

      print("# Remove rpl,rps,point genes")
      rpl_genes = adata.var_names.str.startswith('Rpl')
      adata=adata[:,~rpl_genes]

      rps_genes = adata.var_names.str.startswith('Rps')
      adata=adata[:,~rps_genes]

      point_genes=np.array([True if "." in var else False for var in adata.var_names])
      adata=adata[:,~point_genes]

      number_genes=np.array([True if len(re.findall('^[0-9]+', var))>0 else False for var in adata.var_names])
      adata=adata[:,~number_genes]

      aggregation=os.path.join(path,"outs/aggregation.csv")
      aggr=pd.read_csv(aggregation,sep=",")
      aggr["idents"]=list(range(1,aggr.shape[0]+1))
      status_dict={str(i):s for i,s in zip(aggr["idents"].values.tolist(),aggr["status"].values.tolist())}

      status=[status_dict[i] for i in  adata.obs["idents"]]
      adata.obs["status"]=status

      adata.raw=adata
      return adata
######################################


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=str,default="output")
    parser.add_argument("--path",type=str,default=None,help="cellranger count and aggr output ")

    parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset clusters in column for downstream")
    parser.add_argument("--column",type=str,default=None,help="the column in metadata using for subset,eg:orig_cluster")

    parser.add_argument("--reverse", action='store_true', default=False,help="reverse for subset if true")
    parser.add_argument("--no_doublet", action='store_true', default=False,help="wether remove doublet")
    args=parser.parse_args()

    if not os.path.exists(args.outdir):
       os.makedirs(args.outdir)
    
    sc.settings.autoshow=False
    sc.settings.figdir=args.outdir
    adata=load_data(args.path,args.column,args.subset,args.reverse)

    adata=adata[adata.obs.n_genes_by_counts > 250,:]
    adata=adata[adata.obs.n_genes_by_counts < 2500,:]
    adata=adata[adata.obs.pct_counts_mt < 6,:]

    print("### Detect Doublet")
    model=DoubletModule(adata,outdir=args.outdir)
    model.detect()
    model.plot_histogram()
    model.set_embedding()
    model.add_embedding()
    adata=model.adata
    if args.no_doublet:
        adata=adata[~adata.obs.predicted_doublets,:]

    print("### save")
    adata.write(os.path.join(args.outdir,'adata.h5ad'), compression='gzip')

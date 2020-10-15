import os
import scanpy as sc
import scirpy as ir
import argparse

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from collections import Counter

class Model(object):
   def __init__(self,path,tcr_path=None,n_top_genes=5000,outdir="./scanpy_result"):
      self._path=path
      self._outdir=outdir
      self._tcr_path=tcr_path

      self._n_top_genes=n_top_genes

      self.graphclust_path=os.path.join(self._path,"outs/analysis/clustering/graphclust/clusters.csv")
      self.km_path=os.path.join(self._path,"analysis/clustering/kmeans_10_clusters/clusters.csv")
      self.mtx_path=os.path.join(self._path,"outs/filtered_feature_bc_matrix")
   
   def _load_data(self):
      adata=sc.read_10x_mtx(self.mtx_path,cache=True)
      orig_cluster=pd.read_csv(self.graphclust_path)
      km_cluster=pd.read_csv(self.km_path)

      print("*** Add original message")
      adata.obs["orig_cluster"]=orig_cluster.Cluster.to_list()
      adata.obs["km_cluster"]=km_cluster.Cluster.to_list()

      cells=adata.obs_names.to_list()
      idents=[cell.split("-")[1] for cell in cells] 
      adata.obs["idents"]=idents

      self._n_cells=adata.shape[0]
      self._n_features=adata.shape[1]
      print("origin adata shape is [{},{}]".format(self._n_cells,self._n_features))
      adata=self._remove(adata)  # remove dropout features

      if self._tcr_path is not None:
            if os.path.exists(self._tcr_path):
               adata.uns["tcr_path"]=self._tcr_path
               try:
                  adata_tcr=ir.read_10x_vdj(self._tcr_path)
                  print("*** Merge adata with tcr")
                  ir.pp.merge_with_tcr(adata,adata_tcr)
                  
                  print("*** chain pairing")
                  ir.tl.chain_pairing(adata)

                  ir.pl.group_abundance(adata, groupby="chain_pairing", target_col="orig_cluster")
                  plt.savefig(os.path.join(args.outdir,"chain_pairing.pdf"))
                  plt.close()

                  print("Fraction of cells with more than one pair of TCRs: {:.2f}".format(np.sum(
                     adata.obs["chain_pairing"].isin(["Extra beta", "Extra alpha", "Two full chains"])) / adata.n_obs))
               except:
                  print("Add TCR Invalid!!!")
      
      #if self._subset is not None:
      #  adata=self._select(adata)
      #self.adata=adata
      adata.uns["path"]=self._path
      return adata
   
   def subset_by_cluster(self,adata,subset,col_use="orig_cluster",invert=False):
      metadata=adata.obs
      assert col_use in metadata.columns
      #sets=metadata[col_use].values.tolist()
      adata.obs[col_use].astype("category")
      #sets=[str(s) for s in sets]
      subset=[str(s) for s in subset]
      cols=np.unique(adata.obs[col_use].values.tolist())
      print(subset)
      if invert:
          #subset=[s for s in cols if s not in subset]
          adata_subset=adata[~adata.obs[col_use].isin(subset)]
      else:
          adata_subset=adata[adata.obs[col_use].isin(subset)]
      print("after subset,adata shape is [{},{}]".format(adata_subset.shape[0],adata_subset.shape[1]))
      return adata_subset

   def _remove(self,adata):
      print("remove dropout features")
      mito_genes = adata.var_names.str.startswith('MT-')
      adata=adata[:,~mito_genes]
      rpl_genes = adata.var_names.str.startswith('RPL')
      adata=adata[:,~rpl_genes]

      rps_genes = adata.var_names.str.startswith('RPS')
      adata=adata[:,~rps_genes]
      return adata
   
   def _process(self,adata,batch_key=None):
      adata.obs['n_counts'] =adata.X.sum(axis=1).A1
      print("**** data transform ****")
      sc.pp.normalize_total(adata, target_sum=1e4)
      sc.pp.log1p(adata)
      adata.raw = adata

      sc.pp.highly_variable_genes(adata,n_top_genes=self._n_top_genes)
      #sc.pp.regress_out(adata, ['n_counts'])
      adata = adata[:, adata.var.highly_variable]
      sc.pp.regress_out(adata, ['n_counts'])
      print("**** scale data ****")
      sc.pp.scale(adata, max_value=10)
      
      idx=adata.var["highly_variable"]
      var_genes=adata.var_names[idx]
      var_genes=var_genes.to_list()
      adata.uns["HVG"]=var_genes
      if batch_key is not None:
          assert batch_key in adata.obs.columns
          print("**** correct batch effect ****")
          sc.pp.combat(adata,key=batch_key)
          #sc.external.pp.mnn_correct(adata,var_subset=var_genes,batch_key=key)
          #sc.external.pp.bbknn(adata,batch_key=key,approx=True)
      return adata
   
   def _reduction(self,adata,resolution=1.2):
      print("*** Principal component analysis ***")
      print("**** run pca ****")
      sc.tl.pca(adata, svd_solver='arpack')
      sc.pl.pca_variance_ratio(adata, log=True,save=True)
      print("**** run neighborhood graph ****")
      sc.pp.neighbors(adata, n_neighbors=20, n_pcs=30)

      print("**** run tsne ****")
      sc.tl.tsne(adata,n_pcs=30,learning_rate=200)

      print("**** clustering ****")
      sc.tl.leiden(adata,resolution=resolution)
      sc.tl.louvain(adata,resolution=resolution)
      
      print("**** run umap ****")
      sc.tl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
      sc.pl.paga(adata, plot=False)
      sc.tl.umap(adata, init_pos='paga')

      return adata
   
   def _preload(self,pre_path=None):
      adata=sc.read_h5ad(pre_path)
      return adata
   
   def _FindMarkers(self,adata,groupby="leiden"):
      sc.tl.rank_genes_groups(adata,groupby, method='t-test',n_genes=200,corr_method="bonferroni")
      return adata

   def __str__(self):
      output="The class :  scanpy object"+"\n"
      output=output+"The matrix path : {}".format(self.mtx_path) + "\n"
      output=output+"The graph cluster path : {}".format(self.graphclust_path) + "\n"
      output=output+"Select : {} top genes".format(self._n_top_genes) + "\n"
      if self._tcr_path is not None:
         output=output+"The graph TCR path : {}".format(self._tcr_path) + "\n"

      return output


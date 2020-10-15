import os
import numpy as np
import pandas as pd
import scipy
import scvi
import re
from scipy import stats
from typing import Tuple

import os
import scanpy as sc
import scirpy as ir
import scanpy.external as sce
import argparse

import seaborn as sns
import matplotlib.pyplot as plt
import harmonypy as hm
from collections import Counter

from scvi.dataset import AnnDatasetFromAnnData
from scvi.inference import UnsupervisedTrainer
from scvi.models.vae import VAE

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

def sample_adata(adata,ratios=[0.5],column="idents"):
    metadata=adata.obs
    assert column in metadata.columns
    n_column=len(metadata[column].unique())
    column_dict=dict(Counter(metadata[column]))
    if len(ratios)==1:
        ratios=ratios*n_column
    cells=[np.random.choice(metadata.cells[metadata[column].isin([str(c)])].to_list(),size=int(int(column_dict[str(c)])*ratios[i]),replace=False) \
            for i,c in  enumerate(metadata[column].unique())]
    print(len(cells))
    barcodes=[]
    for cell in cells:
        barcodes.extend(cell)
    print("# the size of barcodes : {}".format(len(barcodes)))
    adata=subset_by_cell_feature(adata,cells=barcodes)
    return adata

def batch_correct(adata,batch_key="idents",method="combat"):
    print("# correct batch effect")
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

def process(adata,key=None,method="combat",use_hv=False):
      # method : harmony,combat
      print("# data transform")
      sc.pp.normalize_total(adata, target_sum=1e4)
      sc.pp.log1p(adata)

      sc.pp.highly_variable_genes(adata,n_top_genes=2000)
      #adata.raw = adata
      if use_hv:
          adata = adata[:, adata.var.highly_variable]
      print("# scale data")
      sc.pp.scale(adata, max_value=10)

      if key is not None:
          assert key in adata.obs.columns
          adata=batch_correct(adata,batch_key=key,method=method)
      return adata

def reduction(adata,resolution=1.5,use_rep="X_hamony",batch_correct=False):

      # use_rep : X_pca,X_harmony,X_scvi...
      #if use_rep=="X_pca":
      #    if use_rep not in adata.obsm_keys():
      #        sc.tl.pca(adata,n_comps=50,svd_solver='arpack')
      if not batch_correct:
          sc.tl.pca(adata,n_comps=50,svd_solver='arpack')
      assert use_rep in adata.obsm_keys()
      print("# Run neighborhood graph ")
      sc.pp.neighbors(adata, n_neighbors=20,use_rep=use_rep)

      print("# Run tsne")
      sc.tl.tsne(adata,use_rep=use_rep,learning_rate=200)

      print("# Clustering")
      sc.tl.leiden(adata,resolution=resolution)
      sc.tl.louvain(adata,resolution=resolution)

      print("# Run umap")
      sc.tl.paga(adata,groups="leiden", use_rna_velocity=False)  # remove `plot=False` if you want to see the coarse-grained graph
      sc.pl.paga(adata, plot=False)
      sc.tl.umap(adata, init_pos='paga')

      return adata


def load_data(path,column=None,subset=None,reverse=False,sample_column="idents"):
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
      #adata.raw=adata
      return adata

def preprocess(adata):
        #adata.obs['n_counts'] =adata.X.sum(axis=1).A1
        #adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1
        #adata.var['mt'] = adata.var_names.str.startswith('MT-')
        #adata.var['rpl'] = adata.var_names.str.startswith('RPL')
        #adata.var['rps'] = adata.var_names.str.startswith('RPS')
        #sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)

        raw=adata.copy()
        print("# Remove useless genes")
        mito_genes = adata.var_names.str.startswith('mt-')
        adata=adata[:,~mito_genes]

        rpl_genes = adata.var_names.str.startswith('Rpl')
        adata=adata[:,~rpl_genes]

        rps_genes = adata.var_names.str.startswith('Rps')
        adata=adata[:,~rps_genes]

        point_genes=np.array([True if "." in var else False for var in adata.var_names])
        adata=adata[:,~point_genes]

        number_genes=np.array([True if len(re.findall('^[0-9]+', var))>0 else False for var in adata.var_names])  
        adata=adata[:,~number_genes]

        print("# After remove useless genes, adata shape is [{},{}]".format(adata.shape[0],adata.shape[1]))
        adata_original = adata.copy()
        print("# Data transform ")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        adata.raw = raw[:,adata.var_names]
        sc.pp.highly_variable_genes(adata,n_top_genes=2000)
        highly_variable_genes = adata.var["highly_variable"]
        #adata = adata[:, highly_variable_genes]

        print("# Scale data")
        sc.pp.regress_out(adata, ["n_counts", "pct_counts_mt"])
        sc.pp.scale(adata, max_value=10)
        adata_original = adata_original[:, highly_variable_genes]
        print("{} highly variable genes".format(highly_variable_genes.sum()))
        
        # adata_original only for scvi
        return adata,adata_original


def load_tcr(adata,tcr_path):
        adata.uns["tcr_path"]=tcr_path
        try:
            adata_tcr=ir.read_10x_vdj(tcr_path)
            print("# Merge adata with tcr")
            ir.pp.merge_with_tcr(adata,adata_tcr)

            print("# Chain pairing")
            ir.tl.chain_pairing(adata)
            print("Fraction of cells with more than one pair of TCRs: {:.2f}".format(np.sum(
                     adata.obs["chain_pairing"].isin(["Extra beta", "Extra alpha", "Two full chains"])) / adata.n_obs))
        except:
            print("Add TCR Invalid!!!")

        return adata

def compute_scvi_latent(
    adata: sc.AnnData,
    n_latent: int = 10,
    n_epochs: int = 100,
    lr: float = 1e-3,
    use_batches: bool = False,
    use_cuda: bool = False,
) -> Tuple[scvi.inference.Posterior, np.ndarray]:
    """Train and return a scVI model and sample a latent space

    :param adata: sc.AnnData object non-normalized
    :param n_latent: dimension of the latent space
    :param n_epochs: number of training epochs
    :param lr: learning rate
    :param use_batches
    :param use_cuda
    :return: (scvi.Posterior, latent_space)
    """
    # Convert easily to scvi dataset
    scviDataset = AnnDatasetFromAnnData(adata)

    # Train a model
    vae = VAE(
        scviDataset.nb_genes,
        n_batch=scviDataset.n_batches * use_batches,
        n_latent=n_latent,
    )
    trainer = UnsupervisedTrainer(vae, scviDataset, train_size=1.0, use_cuda=use_cuda)
    trainer.train(n_epochs=n_epochs, lr=lr)
    ####

    # Extract latent space
    posterior = trainer.create_posterior(
        trainer.model, scviDataset, indices=np.arange(len(scviDataset))
    ).sequential()

    latent, _, _ = posterior.get_latent()

    return posterior, latent

def runscVi(adata,
        n_hidden=128,
        n_layers=1,
        n_epochs=400,
        n_latent=10,
        lr=1e-3,
        batch_key="idents",
        use_cuda=False):
    print("Training with scvi...")
    sce.pp.scvi(adata,n_hidden=n_latent,
        n_latent=n_latent,
        n_layers=n_layers,
        n_epochs=n_epochs,
        lr=lr,
        batch_key=batch_key,
        use_cuda=use_cuda)
    return adata

def scviUpdate(adata,adata_original):
    adata.obsm['X_scvi']=adata_original.obsm['X_scvi']
    adata.obsm['X_scvi_denoised']=adata_original.obsm['X_scvi_denoised']
    adata.obsm['X_scvi_sample_rate']=adata_original.obsm['X_scvi_sample_rate']
    return adata

def run_harmony(adata,batch_key="status",max_iter=10):
    data=adata.copy()
    print("# Run PCA")
    sc.tl.pca(data,n_comps=50,svd_solver='arpack')
    X=adata.X
    print("# The X matrix for harmony shape is : [{},{}]".format(X.shape[0],X.shape[1]))
    meta_data=data.obs
    assert batch_key in meta_data.columns
    vars_use=[batch_key]
    print("# Run harmony")
    ho = hm.run_harmony(X, meta_data, vars_use,max_iter_kmeans=max_iter)
    X_harmony=ho.Z_corr.transpose()

    print("# Add harmony ***")
    data.uns["harmony"]={}
    data.uns["harmony"]["params"]=ho.__dict__
    data.obsm["X_harmony"]=X_harmony
    return data


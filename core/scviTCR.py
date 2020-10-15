import os
import scanpy as sc
import scirpy as ir
import scanpy.external as sce
import argparse

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt

from collections import Counter

import scvi
from scvi.dataset import AnnDatasetFromAnnData
from scvi.inference import UnsupervisedTrainer
from scvi.models.vae import VAE

from .tl import run_harmony
from .scrubletDoublet import DoubletModule
from .utils import subset_by_column,subset_by_cell_feature

sc.settings.autoshow=False
class scviTCR:
    def __init__(self,path,tcr_path=None,batch_key="idents"):
        self.path=path
        self.tcr_path=tcr_path
        self.batch_key=batch_key
        
        self.graphclust_path=os.path.join(self.path,"outs/analysis/clustering/graphclust/clusters.csv")
        self.km_path=os.path.join(self.path,"outs/analysis/clustering/kmeans_10_clusters/clusters.csv")
        self.mtx_path=os.path.join(self.path,"outs/filtered_feature_bc_matrix")

    def _load_data(self):
        print("# Loading 10X data from : {}".format(self.mtx_path))
        adata=sc.read_10x_mtx(self.mtx_path,cache=True)
        adata.var_names_make_unique()
        orig_cluster=pd.read_csv(self.graphclust_path)
        km_cluster=pd.read_csv(self.km_path)

        print("# Add original message")
        print("# Loading graph base clusters from : {}".format(self.graphclust_path))
        print("# Loading kmeans base clusters from : {}".format(self.km_path))
        adata.obs["orig_cluster"]=[str(i) for i in orig_cluster.Cluster.to_list()]
        adata.obs["km_cluster"]=[str(i) for i in km_cluster.Cluster.to_list()]

        cells=adata.obs_names.to_list()
        adata.obs["cells"]=cells
        idents=[cell.split("-")[1] for cell in cells]
        adata.obs["idents"]=idents

        self.n_cells=adata.shape[0]
        self.n_features=adata.shape[1]
        print("# Origin adata shape is [{},{}]".format(self.n_cells,self.n_features))
        print("# Preprocess ")
        adata,adata_original=self._preprocess(adata)
        if self.tcr_path is not None and os.path.exists(self.tcr_path):
            adata=self._load_tcr(adata,self.tcr_path)
        return adata,adata_original



    def _preprocess(self,adata):
        adata.obs['n_counts'] =adata.X.sum(axis=1).A1
        adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        adata.var['rpl'] = adata.var_names.str.startswith('RPL')
        adata.var['rps'] = adata.var_names.str.startswith('RPS')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)

        print("# Remove useless genes")
        mito_genes = adata.var_names.str.startswith('MT-')
        adata=adata[:,~mito_genes]
        
        rpl_genes = adata.var_names.str.startswith('RPL')
        adata=adata[:,~rpl_genes]

        rps_genes = adata.var_names.str.startswith('RPS')
        adata=adata[:,~rps_genes]

        point_genes=np.array([True if "." in var else False for var in adata.var_names])
        adata=adata[:,~point_genes]

        adata_original = adata.copy()
        print("# Data transform ****")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata,n_top_genes=5000)
        highly_variable_genes = adata.var["highly_variable"]
        adata = adata[:, highly_variable_genes]

        adata.raw = adata
        print("# Scale data")
        sc.pp.regress_out(adata, ["n_counts", "mt"])
        sc.pp.scale(adata, max_value=10)
        adata_original = adata_original[:, highly_variable_genes]
        print("{} highly variable genes".format(highly_variable_genes.sum()))

        return adata,adata_original

    def _load_tcr(self,adata,tcr_path):
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

    def _runscVI(self):
        adata,adata_original=self._load_data()
        adata_original=runscVi(adata_adata_original)
        print("# Update adata obsm")
        adata.obsm['X_scvi']=adata_original.obsm['X_scvi']
        adata.obsm['X_scvi_denoised']=adata_original.obsm['X_scvi_denoised']
        adata.obsm['X_scvi_sample_rate']=adata_original.obsm['X_scvi_sample_rate']
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
        n_epochs=500,
        n_latent=10,
        lr=1e-3,
        batch_key="idents"):
    print("Training with scvi...")
    sce.pp.scvi(adata,n_hidden=n_latent,
        n_latent=n_latent,
        n_layers=n_layers,
        n_epochs=n_epochs,
        lr=lr,
        batch_key=batch_key,
        use_cuda=False)
    return adata





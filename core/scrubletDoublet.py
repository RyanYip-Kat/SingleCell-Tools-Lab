import os
import scanpy as sc
import numpy as np
import pandas as pd
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




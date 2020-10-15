import os
import scanpy as sc
import doubletdetection
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt

class DoubletDT():
    def __init__(self,outdir,n_iters=50,save=True):
        self.outdir=outdir
        self.n_iters=n_iters
        self.save=save
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    def detect(self,path,copy=True):
        try:
            adata=sc.read_h5ad(path)
        except:
            adata=sc.read_10x_mtx(path)

        raw_counts=adata.raw.X
        zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
        print("There are {} gene have zero expression".format(sum(zero_genes)))
        raw_counts = raw_counts[:, ~zero_genes]

        print("*** Create doubletdetection BoostClassifier ***")
        clf = doubletdetection.BoostClassifier(n_iters=self.n_iters, 
                use_phenograph=False, standard_scaling=True)

        model=clf.fit(raw_counts)
        doublets=model.predict(p_thresh=1e-16, voter_thresh=0.5)
        doublets=[str(int(i)) for i in doublets]
        ds={"0":"Single","1":"doublet"}
        doublets_predict=[ds[i] for i in doublets]
        print("*** Add predicted doublets on scanpy adata object ***")
        
        print("*** Plot metrics ***")
        doubletdetection.plot.convergence(model, save=os.path.join(self.outdir,'convergence_test.pdf'), show=False, p_thresh=1e-16, voter_thresh=0.5)
        doubletdetection.plot.threshold(model, save=os.path.join(self.outdir,'threshold_test.pdf'), show=False, p_step=6)
        ###################
        if copy:
            data=adata.copy()
        else:
            data=adata
        data.obs["doublets"]=doublets_predict

        if self.save:
            data.write(os.path.join(self.outdir,"adata_doubletdetection.h5ad"))
        return data







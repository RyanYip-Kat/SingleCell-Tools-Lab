import episcanpy as epi
import scanpy.api as sc
import os

import numpy as np
from scipy.io import mmread
import pandas as pd
import anndata as ad


def read_ATAC(path):
    # read the mtx file, the tsv file cell_names and the bed file for var_names"

    files=os.listdir(path)
    assert "matrix.mtx" in files and "barcodes.tsv" in files and "peaks.bed" in files
    matrix_file=os.path.join(path,"matrix.mtx")
    cell_file=os.path.join(path,"barcodes.tsv")
    var_file=os.path.join(path,"peaks.bed")
    mat = mmread(matrix_file)
    mat = mat.toarray()
    mat = np.matrix(mat.transpose())
    
    with open(cell_file) as f:
        barcodes = f.readlines()
        barcodes = [x[:-1] for x in barcodes]
        
    with open(var_file) as f:
        var_names = f.readlines()
        var_names = ["_".join(x[:-1].split('\t')) for x in var_names]

    adata = ad.AnnData(mat, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=var_names))
    adata.uns['omic'] = 'ATAC'

    return(adata)



class  epiScpyQC:
    def __init__(self,path,qc_metrix_path):
        self._path=path
        self._qc_metrix_path=qc_metrix_path
        if not os.path.exists(self._qc_metrix_path):
            os.makedirs(self._qc_metrix_path)

        self.adata=read_ATAC(path)

    def cal_QC(self):
       # preliminary filtering to remove cell barcode containing no open features
       epi.pp.filter_cells(self.adata, min_features=1)
       # preliminary filtering to remove features that are not sequenced in any cell
       epi.pp.filter_features(self.adata, min_cells=1)

       print("### plot cell coverage ###")
       epi.pp.coverage_cells(self.adata, binary=True, log=False, bins=50,
               threshold=1000, save=os.path.join(self._qc_metrix_path,'coverage_cells.png'))
       
       epi.pp.coverage_cells(self.adata, binary=True, log=10, bins=50,
               threshold=1000, save=os.path.join(self._qc_metrix_path,'coverage_cells_log10.png'))

       print("### plot feature coverage ###")
       # plot feature coverage
       epi.pp.coverage_features(self.adata, binary=True, log=False, bins=50,
               threshold=40, save=os.path.join(self._qc_metrix_path,'coverage_cells.png'))
       epi.pp.coverage_features(self.adata, binary=True, log=10, bins=50,
               threshold=40, save=os.path.join(self._qc_metrix_path,'coverage_cells_log10.png'))

       epi.pp.density_features(self.adata, hist=True,save=os.path.join(self._qc_metrix_path,"density_features.png"))
       epi.pp.cal_var(self.adata,show=False,save=os.path.join(self._qc_metrix_path,"cal_var.png"))


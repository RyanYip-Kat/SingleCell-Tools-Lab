import scanpy as sc
import numpy as np
import pandas as pd
import argparse
import os
######################################


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=str,default="output")
    parser.add_argument("--data",type=str,default=None,help="cellranger count and aggr output ")

    parser.add_argument("--min_genes",type=int,default=2000,help="min genes for filtering cells")
    parser.add_argument("--max_genes",type=int,default=6000,help="max genes for filtering cells")

    parser.add_argument("--pct_mt", type=float, default=6.0,help="the mt gene pecent")
    args=parser.parse_args()

    if not os.path.exists(args.outdir):
       os.makedirs(args.outdir)
    
    sc.settings.autoshow=False
    sc.settings.figdir=args.outdir
    adata=sc.read_h5ad(args.data)
    print("Before QC,the data shape is {},{}".format(adata.shape[0],adata.shape[1]))
    sc.pp.filter_cells(adata,min_genes=args.min_genes)
    sc.pp.filter_cells(adata,max_genes=args.max_genes)
    adata=adata[adata.obs.pct_counts_mt < args.pct_mt,:]
    print("After QC,the data shape is {},{}".format(adata.shape[0],adata.shape[1]))
    
    metadata=adata.obs.copy()
    metadata["barcode"]=metadata.index
    meta=metadata[["barcode","Donor","Sample","Status"]]
    #meta=metadata[["barcode","n_genes_by_counts"]]
    meta.to_csv(os.path.join(args.outdir,"metadata.csv"),sep=",",index=False)


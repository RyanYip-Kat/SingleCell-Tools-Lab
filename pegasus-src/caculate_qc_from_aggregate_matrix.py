import scanpy as sc
import pegasus as pg
import numpy as np
import argparse
import os


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--h5sc",type=str,default=None,help="cellranger count and aggr output ")
    parser.add_argument("--outdir",type=str,default="output")

    args=parser.parse_args()

    if not os.path.exists(args.outdir):
       os.makedirs(args.outdir)

    sc.settings.autoshow=False
    sc.settings.figdir=args.outdir

    adata=pg.read_input(args.h5sc)
    adata.obs['n_counts'] =adata.X.sum(axis=1).A1
    adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['rpl'] = adata.var_names.str.startswith('Rpl')
    adata.var['rps'] = adata.var_names.str.startswith('Rps')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
            jitter=True,multi_panel=True,size=0.25,show=False,save="_QC.pdf",figsize=[16,10])

    filename=os.path.dirname(args.h5sc)
    adata.write(os.path.join(filename,"adata.h5ad"))

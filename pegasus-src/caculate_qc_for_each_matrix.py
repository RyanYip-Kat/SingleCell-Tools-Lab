import scanpy as sc
import numpy as np
import pandas as pd
import argparse
import os


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--meta",type=str,default=None,help="cellranger count and aggr output ")
    parser.add_argument("--outdir",type=str,default="output")

    args=parser.parse_args()

    if not os.path.exists(args.outdir):
       os.makedirs(args.outdir)

    sc.settings.autoshow=False
    sc.settings.figdir=args.outdir

    meta=pd.read_csv(args.meta)
    #files=meta.Location.values.tolist()
    #samples=meta.Sample.values.tolist()
    #i=1
    for i in range(meta.shape[0]):
        df=meta.loc[i,:]
        sample=df.Sample
        h5=df.Location
        print("### Loading {} file from :  {}".format(sample,h5))
        adata=sc.read_10x_h5(h5)
        print("### The size of {}  is {},{}".format(sample,adata.shape[0],adata.shape[1]))
        adata.obs['n_counts'] =adata.X.sum(axis=1).A1
        adata.obs["n_genes"]=np.sum(adata.X>0,axis=1).A1
        adata.var['mt'] = adata.var_names.str.startswith('mt-')
        adata.var['rpl'] = adata.var_names.str.startswith('Rpl')
        adata.var['rps'] = adata.var_names.str.startswith('Rps')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','rpl','rps'], percent_top=None, log1p=False, inplace=True)
        sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
              jitter=True,multi_panel=True,size=0.25,show=False,save="_"+sample+"_QC.pdf",figsize=[16,10])

        filename=os.path.join(args.outdir,sample)
        if not os.path.exists(filename):
            os.makedirs(filename)
        barcode=[cell.split("-")[0]+"-"+str(df.Donor) for cell in adata.obs_names]
        pbarcode=[str(df.Sample)+"-"+cell.split("-")[0] for cell in adata.obs_names]
        adata.obs_names=barcode
        adata.obs["Sample"]=df.Sample
        adata.obs["Status"]=df.Status
        adata.obs["Donor"]=df.Donor
        adata.obs["barcode"]=barcode
        adata.obs["pbarcode"]=pbarcode
        adata.write(os.path.join(filename,"adata.h5ad"))
        

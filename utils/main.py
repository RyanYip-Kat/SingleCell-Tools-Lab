import argparse
import os
import scanpy as sc
import numpy as np
import sys
import pandas as pd

sys.path.append("../")

from core.functions import reduction,process
from core.functions import subset_by_column,subset_by_cell_feature,load_data,load_tcr,runscVi,scviUpdate,preprocess,batch_correct,run_harmony
from core.utils import AddMetaData
from core.scrubletDoublet import DoubletModule
from collections import Counter

if __name__=="__main__":
    print("")
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=str,default="output")
    parser.add_argument("--path",type=str,default=None,help="10x cellranger count or aggr,reanalysis output")
    parser.add_argument("--subset",nargs="+",type=str,default=None,help="subset clusters in column for downstream")
    parser.add_argument("--column",type=str,default=None,help="the column in metadata using for subset,eg:orig_cluster")
    parser.add_argument("--reverse", action='store_true', default=False,help="reverse for subset if true")
    parser.add_argument("--batch_key",type=str,default="idents")
    parser.add_argument("--n_jobs",type=int,default=12)
    
    parser.add_argument("--no_doublet", action='store_true', default=False,help="whether remove doublet cells")
    parser.add_argument("--scvi", action='store_true', default=False,help="whether to run scvi")
    parser.add_argument("--gpu", action='store_true', default=False,help="whether to use gpu for scvi")

    parser.add_argument("--batch_correct", action='store_true', default=False,help="batch correct")
    parser.add_argument("--batch_method",type=str,default="combat",choices=["combat","mnn","bbknn","harmony"])
    parser.add_argument("--tcr_path",type=str,default=None)
    parser.add_argument("--use_rep",type=str,default="X_harmony",choices=["X_pca","X_harmony","X_scvi"],help="The embedding for tsne,umap,clusters")

    args=parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    
    print("")
    sc.settings.figdir=args.outdir
    sc.settings.autoshow=False
    sc.settings.autosave=True
    sc.settings.cachedir='./cache'

    sc.settings.cache_compression="gzip"
    sc.settings.n_jobs=args.n_jobs
    sc.settings.file_format_figs="pdf"
    
    use_rep=args.use_rep
    aggregation=os.path.join(args.path,"outs/aggregation.csv")
    aggr=pd.read_csv(aggregation,sep=",")
    aggr["idents"]=list(range(1,aggr.shape[0]+1))
    status_dict={str(i):s for i,s in zip(aggr["idents"].values.tolist(),aggr["status"].values.tolist())}

    print("# Laoding data from path : {}".format(args.path))
    adata=load_data(args.path,args.column,args.subset,args.reverse)
    if args.tcr_path is not None:
        adata=load_tcr(adata,args.tcr_path)
    
    #adata=adata[adata.obs.n_genes_by_counts > 200,:]
    #adata=adata[adata.obs.n_genes_by_counts < 2500,:]
    #adata=adata[adata.obs.pct_counts_mt < 6,:]
    #print(adata.shape)

    status=[status_dict[i] for i in  adata.obs["idents"]]
    adata.obs["status"]=status
    Counter(status)

    #df=pd.read_csv("/home/ye/Work/Python/SingleCell/Project/dr2020/20200727/RBC.csv")
    #cells=df.iloc[:,0].values.tolist()
    #adata=subset_by_cell_feature(adata,cells=cells,invert=True)

    print("# Preprocess adata")
    adata,adata_original=preprocess(adata)
    filename=os.path.join(args.outdir,"adata_raw.h5ad")
    adata_original.write(filename)

    print("### Detect Doublet")
    model=DoubletModule(adata,outdir=args.outdir)
    model.detect()
    model.plot_histogram()
    model.set_embedding()
    model.add_embedding()
    adata=model.adata
    if args.no_doublet:
        adata=adata[~adata.obs.predicted_doublets,:]
        #adata_original=adata_original[~adata.obs.predicted_doublets,:]

    if args.batch_correct:
        print("# Use_rep is {},and the batch correct method is {}".format(use_rep,args.batch_correct))
        adata=batch_correct(adata,batch_key=args.batch_key,method=args.batch_method)

    #if use_rep=="X_harmony" and args.batch_key is not None:
    #    print("# Run harmony")
    #    adata=run_harmony(adata,batch_key=args.batch_key)

    if args.scvi:
        print("# Run scvi")
        cells=adata.obs_names
        adata_original=subset_by_cell_feature(adata_original,cells=cells)
        adata_original=runscVi(adata_original,batch_key=args.batch_key,use_cuda=args.gpu)
   
        filename=os.path.join(args.outdir,"adata_scvi.h5ad")
        print("# Save original data into : {}".format(filename))
        adata_original.write(filename)

        adata=scviUpdate(adata,adata_original)
        use_rep="X_scvi"
    
    print("# Reductions")
    adata=reduction(adata,use_rep=use_rep,resolution=1.2)
    
    print("# FindMarkers")
    sc.tl.rank_genes_groups(adata,groupby="leiden", method='t-test',n_genes=200,corr_method="bonferroni")
    filename=os.path.join(args.outdir,"adata.h5ad")
    print("# Save adata into : {}".format(filename))
    adata.write(filename)
    print("# Done!")




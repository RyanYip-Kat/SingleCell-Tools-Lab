import argparse
import os
import scanpy as sc
import numpy as np
import pandas as pd
import sys

sys.path.append("../")

from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--n_top",type=int,default=None)
parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default="leiden",help="data")
parser.add_argument("--column2",type=str,default=None,help="data")
parser.add_argument("--invert", action='store_true', default=False)
parser.add_argument("--sample_n",type=int,default=None)
parser.add_argument("--sample_col",type=str,default="idents")


args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
print("*** Loading data")
data=sc.read_h5ad(args.data)
print("# Origin data shape is :[{},{}]".format(data.shape[0],data.shape[1]))

point_genes=np.array([True if "." in var else False for var in data.var_names])
data=data[:,~point_genes]
print("# After remove genes, data shape is :[{},{}]".format(data.shape[0],data.shape[1]))

if args.subset1 is not None:
    data=subset_by_column(data,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        data=subset_by_column(data,args.subset2,args.column2)
adata=data.copy()
metadata=adata.obs
cells=adata.obs_names
metadata["cells"]=cells

if args.sample_n is not None:
    assert args.sample_col is not  None or args.sample_col in metadata.columns
    
    sample_n=args.sample_n
    d_count=Counter(metadata[args.sample_col])
    d_min=min(list(d_count.values()))
    if args.sample_n > d_min:
        sample_n=d_min
    
    cells=[np.random.choice(metadata.cells[metadata[args.sample_col].isin([str(i)])].to_list(),size=sample_n,replace=False) \
            for i in  metadata[args.sample_col].unique()]
    print(len(cells))
    barcodes=[]
    for cell in cells:
        barcodes.extend(cell)
    print("# the size of barcodes : {}".format(len(barcodes))) 
    adata=subset_by_cell_feature(adata,cells=barcodes)

    if adata.raw is not None:
        adata=adata.raw.to_adata()

    point_genes=np.array([True if "." in var else False for var in adata.var_names])
    adata=adata[:,~point_genes]
    print("# After remove genes, data shape is :[{},{}]".format(adata.shape[0],adata.shape[1]))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    if args.n_top is not None:
        sc.pp.highly_variable_genes(adata,n_top_genes=args.n_top)
        adata = adata[:, adata.var.highly_variable]

else:
    if adata.raw is not None:
        adata=adata.raw.to_adata()

    point_genes=np.array([True if "." in var else False for var in adata.var_names])
    adata=adata[:,~point_genes]
    print("# After remove genes, data shape is :[{},{}]".format(adata.shape[0],adata.shape[1]))

metadata=adata.obs
#cells=[cell.replace("-","") for cell in adata.obs_names]
cells=["SRR"+str(i) for i in range(len(adata.obs_names))]
meta=metadata.copy()
meta["cell_id"]=cells
#meta=meta[["status","idents","cell_id","label_fine"]]


filename=os.path.join(args.outdir,"cell_metadata.csv")
meta.to_csv(filename,sep=",",index=False)

filename=os.path.join(args.outdir,"gene_metadata.csv")
gene_meta=pd.DataFrame(adata.var_names,columns=["symbol"])
gene_meta.to_csv(filename,sep=",",index=False)

X=adata.X.toarray()
print(X.shape)
#cells=[cell.replace("-","") for cell in adata.obs_names]
matrix=pd.DataFrame(X.transpose(),columns=cells,index=adata.var_names)

filename=os.path.join(args.outdir,"linear_gene_expression_matrix.tsv")
print("*** Expression matrix output in {}".format(filename))
matrix.to_csv(filename,sep="\t")

print("*** Done! ***")

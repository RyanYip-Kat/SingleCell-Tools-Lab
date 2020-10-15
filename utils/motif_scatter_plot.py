import scanpy as sc
import argparse
import sys
import os
import seaborn as sns


sys.path.append("../")

from core.utils import get_rank_group_genes,subset_by_column

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--data",type=str,default=None,help="data")
parser.add_argument("--motifs",nargs="+",type=str,default=None)
parser.add_argument("--subset1",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--subset2",nargs="+",type=str,default=None,help="subset")
parser.add_argument("--column1",type=str,default="leiden",help="data")
parser.add_argument("--column2",type=str,default=None,help="data")
parser.add_argument("--basic",type=str,default="Scanpy_UMAP__(highly_variable_genes)",\
        choices=["SCENIC_AUC_UMAP","SCENIC_AUC_t-SNE","Scanpy_t-SNE_(highly_variable_genes)","aucell"],help="basic")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

adata=sc.read_h5ad(args.data)
if args.subset1 is not None:
    adata=subset_by_column(adata,args.subset1,args.column1)
    if args.subset2 is not None and args.column2 is not None:
        adata=subset_by_column(adata,args.subset2,args.column2)



sc.settings.figdir=args.outdir
sc.set_figure_params(frameon=False, dpi=600, dpi_save=600,figsize=[12,8])

for motif in args.motifs:
    mo='Regulon('+motif+'_(+))'
    print(mo)
    try:
        sc.pl.scatter(adata,basis=args.basic,color=mo,save="_"+motif+".pdf",show=False)
    except:
        print("Please check motif : {} whether in adata's vars or metadata columns".format(motif))


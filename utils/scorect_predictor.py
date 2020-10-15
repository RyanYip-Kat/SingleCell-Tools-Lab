import scanpy as sc
import os
import sys
import argparse
from collections import Counter

sys.path.append("../")
from core.core import Model
from core import scorect_api as ct
from core.ScanpyConfig import *

parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--column",type=str,default="leiden")
parser.add_argument("--reference",type=str,default="hsPBMC_markers.csv")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

sc.settings.figdir=args.outdir
sc.settings.n_jobs=8
print("*** Loading adata from : {}".format(args.data))
adata=sc.read_h5ad(args.data)
#adata.obs[args.column]=adata.obs[args.column].astype("category")


print("*** Loading reference from : {}".format(args.reference))
ref_marker = ct.read_markers_from_file(args.reference)

print("*** rank genes group")
sc.tl.rank_genes_groups(adata,groupby=args.column, method='t-test',n_genes=1000,corr_method="bonferroni")
marker_df = ct.wrangle_ranks_from_anndata(adata)
print(marker_df.head())


# Score cell types for each cluster 
# Let's set parameters first - K represents the number of genes included in the ranking
# m represents the number of bins used to divide the top K genes.
K = 300
m = 5
# Get the background genes - here, all the genes used to run the differential gene expression test
background = adata.raw.var.index.tolist()
# Now run the function
print("*** run celltype score")
ct_pval, ct_score = ct.celltype_scores(nb_bins=m,
                                        ranked_genes=marker_df,
                                        K_top = K,
                                        marker_ref=ref_marker,
                                        background_genes=background)
print("*** run assign celltype")
# Now assign clusters to cell types
cluster_assign = adata.obs[args.column]
celltype_assign = ct.assign_celltypes(cluster_assignment=cluster_assign, ct_pval_df=ct_pval, ct_score_df=ct_score, cutoff=0.1)
# Add to anndata object
adata.obs['scorect'] = celltype_assign
# Let's compare with the true assignment now! 
sc.pl.tsne(adata, color=[args.column,'scorect'], title=['True','Predicted'], cmap='Set2',save="_scorect",show=False)

adata.write(os.path.join(args.outdir,"adata.h5ad"),compression='gzip')






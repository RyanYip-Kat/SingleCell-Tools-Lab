import scanpy as sc
import scanpy.external as sce
import palantir
import matplotlib.pyplot as plt

import numpy as np
import os
import argparse
import pickle

import sys

sys.path.append("../")

from core.utils import subset_by_column,subset_by_cell_feature
from collections import Counter


parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--subset",nargs="+",type=str,default=None)
parser.add_argument("--column",type=str,default="leiden")
parser.add_argument("--genes",type=str,default=None)
parser.add_argument("--cluster",type=str,default="label_fine")
parser.add_argument("--start_cell",type=str,default=None)


args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)


figure=os.path.join(args.outdir,"figure")
if not os.path.exists(figure):
   os.makedirs(figure)


adata=sc.read_h5ad(args.data)
cells=adata.obs_names

if args.column is not None:
    if args.subset is not None:
        adata=subset_by_column(adata,args.subset,args.column)
        sc.pp.normalize_total(adata, target_sum=1e4)

print("# Run palantir in scanpy")
sce.tl.palantir(adata,normalize=False, 
                log_transform=False,
                filter_low=False,
                inplace=True)

#adata.write(os.path.join(args.outdir,"adata_palantir.h5ad"))
pca_projections=adata.uns['palantir_pca_results']['pca_projections']
variance_ratio=adata.uns['palantir_pca_results']['variance_ratio']

#norm_df=adata.uns['palantir_norm_data'] #if normalize :True
ms_data=adata.uns['palantir_ms_data']  # Multi scale data matrix
tsne=adata.uns['palantir_tsne']  #tSNE on diffusion maps


imp_df =adata.uns['palantir_imp_df']  # MAGIC imputation
dm_res=adata.uns['palantir_diff_maps'] # Diffusion components, corresponding eigen values and diffusion operator

#palantir.plot.plot_gene_expression(imp_df, tsne, ['CD34', 'MPO', 'GATA1', 'IRF8'])

print("# Plot diffusion map")
palantir.plot.plot_diffusion_components(tsne, dm_res)  # diffusion map
plt.savefig(os.path.join(figure,"diffusion_map.pdf"))

# plot on clusters
if args.cluster is not None:
    clusters=adata.obs[args.cluster]
    palantir.plot.plot_cell_clusters(tsne, clusters )
    plt.savefig(os.path.join(figure,args.cluster+"_plot.pdf"))


# Palantir can be run by specifying an approxiate early cell. While Palantir automatically determines the terminal states

print("# Run Palantir with start cell")

start_cell=args.start_cell if args.start_cell is not None else np.random.choice(cells)
pr_res = palantir.core.run_palantir(ms_data,start_cell,n_jobs=8)
with open(os.path.join(args.outdir,"pr_res.pkl"),"wb") as f:
  pickle.dump(pr_res,f)


#genes = ['GATA1', 'IRF8']
#gene_trends = palantir.presults.compute_gene_trends( pr_res, imp_df.loc[:, genes])

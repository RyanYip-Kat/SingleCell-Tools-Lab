import argparse
import os
import scanpy as sc
import numpy as np
import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt

from core.tl import reduction,process
from core.core import Model
from core.utils import subset_by_column,subset_by_cell_feature

parser=argparse.ArgumentParser()
parser.add_argument("-o","--outdir",type=str,default=None)
parser.add_argument("-d","--data",type=str,default="/home/ye/Work/BioAligment/10X/output/20200407_aging-XXX")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)


net=Model(args.data)
print(net)

adata=net._load_data()
adata=net._subset(adata,[10,13])

counts_matrix=adata.raw.X

print("### Initialize Scruble")
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)

print("### Detect ,Nomalize,PCA")
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)



print("### Set doublets score")
scrub.call_doublets(threshold=0.25)

print("### Plot doublets score expected and sim")
scrub.plot_histogram()
plt.savefig(os.path.join(args.outdir,"doublets_score_hist.pdf",figsize=(12, 12))) 

print("### Get Embeding")
print("#### Run UMAP")
scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 30, min_dist=0.3))

# # Uncomment to run tSNE - slow
print('#### Running tSNE...')
scrub.set_embedding('tSNE', scr.get_tsne(scrub.manifold_obs_, angle=0.9,verbose=True))

# # Uncomment to run force layout - slow
print('#### Running ForceAtlas2...')
scrub.set_embedding('FA', scr.get_force_layout(scrub.manifold_obs_, n_neighbors=20,n_iter=1000))

print('Done.')

print("Export Embedings")
umap_mat=scrub._embeddings["UMAP"]
tsne_mat=scrub._embeddings["tSNE"]
fa_mat=scrub._embeddings["FA"]

umap_df=pd.DataFrame(umap_mat,index=adata.obs_names,columns=["UMAP_1","UMAP_2"])
tsne_df=pd.DataFrame(tsne_mat,index=adata.obs_names,columns=["tSNE_1","tSNE_2"])
#fa_df=pd.DataFrame(fa_mat,index=adata.obs_names,columns=["FA_1","FA_2"])

umap_df.to_csv(os.path.join(args.outdir,"umap.csv"),sep=",")
tsne_df.to_csv(os.path.join(args.outdir,"tsne.csv"),sep=",")



import argparse
import os
import scanpy as sc
import scanpy.external as sce
import time


parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default=None)
parser.add_argument("--data",type=str,default=None)
parser.add_argument("--embedding",type=str,default="umap")

args=parser.parse_args()
if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("### loading dataset")
adata=sc.read_h5ad(args.data)
t0 = time.time()
sce.exporting.spring_project(adata=adata,
        project_dir=args.outdir,
        embedding_method=args.embedding, 
        subplot_name=args.embedding,
        cell_groupings=["leiden","status","louvain","idents"],
        custom_color_tracks=["leiden","status","louvain","idents"])
print(time.time() - t0)

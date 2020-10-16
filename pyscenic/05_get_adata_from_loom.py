# import dependencies
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import json
import base64
import zlib
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import argparse


parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--loom",type=str,default=None,help="scenic_integrated-output.loom")
parser.add_argument("--sig",type=str,default=None,help="reg.csv")

args=parser.parse_args()

if not os.path.exists(args.outdir):
   os.makedirs(args.outdir)

print("-"*16)
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=150)

print("-"*16)
print("Extract relevant data from the integrated loom file {}".format(args.loom))
# scenic output
lf = lp.connect( args.loom, mode='r', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

# create a dictionary of regulons:
regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
    regulons[i] =  list(r[r==1].index.values)

# cell annotations from the loom column attributes:
cellAnnot = pd.concat(
    [
        pd.DataFrame( lf.ca.CellType, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.ClusterID, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.Leiden_clusters_Scanpy, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
        pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
    ],
    axis=1
)
cellAnnot.columns = [
 'CellType',
 'ClusterID',
 'Leiden_clusters_Scanpy',
 'nGene',
 'nUMI']

# capture embeddings:
dr = [
    pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
]
dr_names = [
    meta['embeddings'][0]['name'].replace(" ","_")
]

# add other embeddings
drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )

for i in range( len(drx.columns) ):
    dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
    dr_names.append( meta['embeddings'][i+1]['name'].replace(" ","_").replace('/','-') )

# rename columns:
for i,x in enumerate( dr ):
    x.columns = ['X','Y']

lf.close()
print("-"*16)

adata = sc.read(args.loom, validate=False)
# drop the embeddings and extra attributes from the obs object
adata.obs.drop( ['Embedding','Embeddings_X','Embeddings_Y','RegulonsAUC'], axis=1, inplace=True )

print("add the embeddings into the adata.obsm object")
for i,x in enumerate( dr ):
    adata.obsm[ 'X_'+dr_names[i] ] = x.as_matrix()

# load the regulons from a file using the load_signatures function
sig = load_signatures(args.sig)
adata = add_scenic_metadata(adata, auc_mtx, sig)

adata.write(os.path.join(args.outdir,"sig_adata.h5ad"))


aucell_adata = sc.AnnData(X=auc_mtx.sort_index())
aucell_adata.obs = adata.obs
aucell_adata.write(os.path.join(args.outdir,"aucell_adata.h5ad"))

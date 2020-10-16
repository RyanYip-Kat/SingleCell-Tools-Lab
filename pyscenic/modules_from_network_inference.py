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
from pyscenic.utils import modules_from_adjacencies
from pyscenic.plotting import plot_binarization
from pyscenic.export import add_scenic_metadata
from pyscenic.cli.utils import load_signatures
import argparse


parser=argparse.ArgumentParser()
parser.add_argument("--outdir",type=str,default="output")
parser.add_argument("--loom",type=str,default=None,help="scenic_integrated-output.loom")
parser.add_argument("--tfs",nargs="+",type=str,default="EBF1",help="")
parser.add_argument("--adj",type=str,default="adj.tsv",help="")

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


adjacencies = pd.read_csv(args.adj, index_col=False, sep='\t')
modules = list(modules_from_adjacencies(adjacencies, exprMat))

print("infer : {} transcription factors".format(len(args.tfs)))
for tf in args.tfs:
    print("pickout modules for : {}".format(tf))
    tf_mods = [ x for x in modules if x.transcription_factor==tf ]
    for i,mod in enumerate( tf_mods ):
        print( f'{tf} module {str(i)}: {len(mod.genes)} genes' )
    
    print( f'{tf} regulon: {len(regulons[tf+"_(+)"])} genes' )
    for i,mod in enumerate( tf_mods ):
        with open( os.path.join(args.outdir,tf+'_module_'+str(i)+'.txt'), 'w') as f:
            for item in mod.genes:
                f.write("%s\n" % item)
    with open( os.path.join(args.outdir,tf+'_regulon.txt'), 'w') as f:
        for item in regulons[tf+'_(+)']:
            f.write("%s\n" % item)



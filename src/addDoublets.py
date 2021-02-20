import argparse
import os
import scanpy as sc
import numpy as np
import sys
import pandas as pd

sys.path.append("../")

from core.utils import AddMetaData
from core.scrubletDoublet import DoubletModule
from collections import Counter

if __name__=="__main__":
    print("")
    parser=argparse.ArgumentParser()
    parser.add_argument("--outdir",type=str,default="output")
    parser.add_argument("--data",type=str,default=None,help="10x cellranger count or aggr,reanalysis output")

    args=parser.parse_args()

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    print("")
    sc.settings.figdir=args.outdir
    sc.settings.autoshow=False
    sc.settings.autosave=True
    sc.settings.cachedir='./cache'

    sc.settings.cache_compression="gzip"
    sc.settings.n_jobs=12
    sc.settings.file_format_figs="pdf"
    print("### Loading data")
    adata=sc.read_h5ad(args.data)
    model=DoubletModule(adata,outdir=args.outdir)
    model.detect()
    model.plot_histogram()
    model.set_embedding()
    model.add_embedding()
    adata=model.adata

    filename=os.path.join(args.outdir,"adata_doublet.h5ad")
    print("# Save adata into : {}".format(filename))
    adata.write(filename)
    print("# Done!")

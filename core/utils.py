import os
import numpy as np
import pandas as pd
import scipy
import scvi
from scipy import stats

import os
import scanpy as sc
import scirpy as ir
import scanpy.external as sce
import argparse

import seaborn as sns
import matplotlib.pyplot as plt
from collections import Counter

from scvi.dataset import AnnDatasetFromAnnData
from scvi.inference import UnsupervisedTrainer
from scvi.models.vae import VAE

def subset_by_column(adata,subset,col_use="orig_cluster",invert=False):
      metadata=adata.obs
      assert col_use in metadata.columns
      #sets=metadata[col_use].values.tolist()
      adata.obs[col_use].astype("category")
      #sets=[str(s) for s in sets]
      subset=[str(s) for s in subset]
      cols=np.unique(adata.obs[col_use].values.tolist())
      print(subset)
      for s in subset:
          if s not in cols:
              raise ValueError("Invalid value in column selected")

      if invert:
          #subset=[s for s in cols if s not in subset]
          adata_subset=adata[~adata.obs[col_use].isin(subset)]
      else:
          adata_subset=adata[adata.obs[col_use].isin(subset)]
      print("after subset,adata shape is [{},{}]".format(adata_subset.shape[0],adata_subset.shape[1]))
      return adata_subset

def subset_by_cell_feature(adata,cells=None,features=None,invert=False):
    data=adata.copy()
    data.obs["cells"]=data.obs_names
    Features=adata.var_names
    print("before subset,there are : {} #cells".format(data.shape[0]))
    if cells is not None:
        if invert:
            data=data[~data.obs["cells"].isin(cells)]
        else:
            data=data[data.obs["cells"].isin(cells)]
    if features is not None:
       if invert:
           feature_idx=[False if feature in features else True for feature in Features]
       else:
           feature_idx=[True if feature in features else False for feature in Features]
       data=data[:,feature_idx]
    print("after subset,there are : {} #cells".format(data.shape[0]))
    return data


def get_rank_group_genes(adata,pval=None,fc=0.25,outdir="./markers",
        n_top=None,back=False,
        key="rank_genes_groups",prefix=""):
    # key : rank_genes_groups or rank_genes_groups_filtered
    assert "rank_genes_groups" in adata.uns_keys()
    groupby=adata.uns[key]["params"]["groupby"]
    G=np.unique(adata.obs[groupby].values)
    table=[]
    for g in G:
        try:
           print("Get markers from : {}".format(g))
           t=sc.get.rank_genes_groups_df(adata,group=g,pval_cutoff=pval,log2fc_min=fc,log2fc_max=10)
           t["cluster"]=g
           t=t.sort_values("logfoldchanges",ascending=False)
           if n_top is not None:
              t=t[:n_top]
           
           t=t[t.logfoldchanges>0]
           t=t[t.pvals<0.05]
           table.append(t)
        except:
            print("group : {} has no rank".format(g))
     
    df=pd.concat(table,axis=0)
    if n_top is not None:
        filename=os.path.join(outdir,prefix+"_"+str(n_top)+"_markers.csv")
    else:
        filename=os.path.join(outdir,prefix+"_markers.csv")
    print("Save markers in {}".format(filename))
    df.to_csv(filename,sep=",",index=False)
    if back:
        return df


def AddMetaData(data,meta,inplace=False):
    if not inplace:
        adata=data.copy()
    else: 
        adata=data
    metadata=adata.obs
    assert meta.shape[0]==metadata.shape[0]
    index=metadata.index.tolist()
    DATA=meta.loc[index,:]  #  loc accept str index
    columns=DATA.columns.tolist()
    for column in columns:
        if column in metadata.columns:
            #new_column="X_"+str(column)
            continue
        else:
            new_column=str(column)
            adata.obs[new_column]=DATA[new_column].values
    if inplace:
        return 
    else:
        return adata




def gene_comparison(adata,gene,by="status"):
    #assert by in ["status","idents"]
    DATA=adata.copy()
    DATA=DATA.raw.to_adata()
    df=DATA.to_df()
    df=df[[gene]]
    df[[by]]=DATA.obs[[by]]
    status=np.unique(df[[by]])

    if len(status)==1:
        raise ValueError("status length must be > = 2")
    elif len(status)==2:
        func=stats.ttest_ind
        x=df[df[by].isin([status[0]])][gene].values
        y=df[df[by].isin([status[1]])][gene].values
        diff,pvals=func(x,y)
        print("t-test,The gene : {}  mean difference between : {} and {} is : {},and the pvals is {}".format(gene,status[0],status[1],diff,pvals))
    elif len(status)>2:
        func=stats.f_oneway
        v_list=[]
        for s in status:
            x=df[df[by].isin([s])][gene].values
            v_list.append(x)

        diff,pvals=func(*v_list)
        print("oneway anova,The gene : {} mean difference across status is :{},and the pvals is {}".format(gene,diff,pvals))

    return diff,pvals


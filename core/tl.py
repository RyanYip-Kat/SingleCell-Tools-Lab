import scanpy as sc
import harmonypy as hm

def process(adata,key=None,method="harmony",use_hv=True):
      # method : harmony,combat
      print("# data transform")
      sc.pp.normalize_total(adata, target_sum=1e4)
      sc.pp.log1p(adata)

      sc.pp.highly_variable_genes(adata,n_top_genes=2000)
      if use_hv:
          adata = adata[:, adata.var.highly_variable]
      print("# scale data")
      sc.pp.scale(adata, max_value=10)

      if key is not None:
          assert key in adata.obs.columns
          print("# correct batch effect with : {}".format(method))
          if method=="harmony":
              adata=run_harmony(adata,batch_key=key)
          elif method=="combat":
              sc.pp.combat(adata,key=key)
          else:
              raise ValueError("Method must be one of hamony or combat!!!")

          #sc.external.pp.mnn_correct(adata,var_subset=var_genes,batch_key=key)
          #sc.external.pp.bbknn(adata,batch_key=key,approx=True)
      return adata

def reduction(adata,resolution=1.5,use_rep="X_hamony"):
      
      # use_rep : X_pca,X_harmony,X_scvi...
      if use_rep=="X_pca":
          if use_rep not in adata.obsm_keys():
              sc.tl.pca(adata,n_comps=50,svd_solver='arpack')
      assert use_rep in adata.obsm_keys()
      print("# Run neighborhood graph ")
      sc.pp.neighbors(adata, n_neighbors=20,use_rep=use_rep)

      print("# Run tsne")
      sc.tl.tsne(adata,use_rep=use_rep,learning_rate=200)

      print("# Clustering")
      sc.tl.leiden(adata,resolution=resolution)
      sc.tl.louvain(adata,resolution=resolution)

      print("# Run umap")
      sc.tl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
      sc.pl.paga(adata, plot=False)
      sc.tl.umap(adata, init_pos='paga')

      return adata



def run_harmony(adata,batch_key="status",max_iter=10):
    data=adata.copy()
    print("# Run PCA")
    sc.tl.pca(data,n_comps=50,svd_solver='arpack')
    #X=data.obsm["X_pca"]
    #if n_pca > X.shape[1]:
    #    n_pca=X.shape[1]
    #X=X[:,:n_pca]
    X=adata.X
    print("# The X matrix for harmony shape is : [{},{}]".format(X.shape[0],X.shape[1]))
    meta_data=data.obs
    assert batch_key in meta_data.columns
    vars_use=[batch_key]
    print("# Run harmony")
    ho = hm.run_harmony(X, meta_data, vars_use,max_iter_kmeans=max_iter)
    X_harmony=ho.Z_corr.transpose()
    
    print("# Add harmony ***")
    data.uns["harmony"]={}
    data.uns["harmony"]["params"]=ho.__dict__
    data.obsm["X_harmony"]=X_harmony
    return data


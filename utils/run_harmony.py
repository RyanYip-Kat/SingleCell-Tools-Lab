import harmonypy as hm

def run_harmony(adata,n_pca=50,batch_key="status",max_iter=10):
    data=data.copy()
    if "X_pca" not in data.obsm_keys():
        print("Run PCA first now!!!")
        sc.tl.pca(data,n_comps=n_pca)
    
    X=data.obsm["X_pca"] 
    if n_pca > X.shape[1]:
        n_pca=X.shape[1]
    X=X[:,:n_pca]
    meta_data=data.obs
    assert batch_key in meta_data.columns
    vars_use=[batch_key]
    ho = hm.run_harmony(X, meta_data, vars_use,max_iter_kmeans=max_iter)
    X_harmony=ho.Z_corr.transpose()
    return X_harmony





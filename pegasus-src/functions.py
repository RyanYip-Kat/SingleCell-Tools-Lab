import numpy as np
from scipy.sparse import issparse

def X_from_rep(data, rep="pca") -> np.array:
    """
    If rep is not mat, first check if X_rep is in data.obsm. If not, raise an error. 
    If rep is None, return data.X as a numpy array
    """
    if rep != "mat":
        rep_key = "X_" + rep
        if rep_key not in data.obsm.keys():
            raise ValueError("Cannot find {0} matrix. Please run {0} first".format(rep))
        return data.obsm[rep_key]
    else:
        return data.X if not issparse(data.X) else data.X.toarray()

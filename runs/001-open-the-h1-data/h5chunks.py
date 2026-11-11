"""Read an .h5ad without loading it: obs columns, var names, and the matrix in row chunks as CSR."""
import h5py
import numpy as np
import pandas as pd
from scipy import sparse


def _read(elem):
    """anndata's own reader handles every on-disk layout (categoricals, nullable strings, old and new versions)."""
    try:
        from anndata.io import read_elem
    except ImportError:  # anndata < 0.11
        from anndata.experimental import read_elem
    return read_elem(elem)


def obs_column(f, name):
    """A missing label (categorical code -1) comes back as NaN and becomes ''. Indexing codes by hand would silently
    return the last category instead, which was run x03's bug."""
    values = pd.Series(_read(f["obs"][name]))
    return values.astype(object).where(values.notna(), "").astype(str).to_numpy()


def obs_names(f):
    return [k for k in f["obs"].keys() if k != f["obs"].attrs.get("_index", "_index")]


def var_names(f):
    var = f["var"]
    key = var.attrs.get("_index", "_index")
    return np.asarray(_read(var[key]), dtype=str)


def shape(f):
    x = f["X"]
    return tuple(int(n) for n in (x.attrs["shape"] if isinstance(x, h5py.Group) else x.shape))


def encoding(f):
    x = f["X"]
    return str(x.attrs.get("encoding-type", "dense" if isinstance(x, h5py.Dataset) else "?"))


def row_chunks(f, size=20000):
    """Yield (start, CSR block) over the rows. Handles CSR and dense storage; CSC is refused rather than read slowly."""
    x = f["X"]
    n, m = shape(f)
    if isinstance(x, h5py.Dataset):
        for s in range(0, n, size):
            yield s, sparse.csr_matrix(x[s:s + size])
        return
    enc = encoding(f)
    if enc != "csr_matrix":
        raise ValueError(f"X is stored as {enc}; row chunking needs CSR")
    indptr = x["indptr"][()]
    for s in range(0, n, size):
        e = min(s + size, n)
        a, b = int(indptr[s]), int(indptr[e])
        yield s, sparse.csr_matrix((x["data"][a:b], x["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, m))

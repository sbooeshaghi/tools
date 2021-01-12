# comparison utils
import os
import pandas as pd
import numpy as np
import scipy
# from openTSNE import TSNE
import openTSNE
from openTSNE.callbacks import ErrorLogger
from sklearn.decomposition import TruncatedSVD
import anndata
from sklearn.metrics.pairwise import manhattan_distances
from sklearn.manifold import TSNE


def import_adata(folder, cr=False):
    '''
        Returns adata with basic structure
    '''
    if cr:
        mat = os.path.join(folder, 'matrix.mtx.gz')
        var = os.path.join(folder, 'barcodes.tsv.gz')
        obs = os.path.join(folder, 'features.tsv.gz')
        bcs = pd.read_csv(var, index_col = 0, header = None, names = ['barcode'])
    else:
        mat = os.path.join(folder, 'genes.mtx')
        var = os.path.join(folder, 'genes.barcodes.txt')
        obs = os.path.join(folder, 'genes.genes.txt')
        bcs = pd.read_csv(var, index_col = 0, header = None, names = ['barcode'])



    if cr:
        gene = pd.read_csv(obs, header = None, index_col = 0, names =['ensembl_id', 'name', 'feature'], sep = '\t')
        gene.index = gene.index.str.slice(0, 18) # slice off the version number from the genes
        adata = anndata.AnnData(X=scipy.io.mmread(mat).tocsr(), obs=gene, var=bcs)
        adata.var.index = adata.var.index.str.slice(0,16,1) # slice the dump -1 off of the barcode
    else:

        gene = pd.read_csv(obs, header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
        gene.index = gene.index.str.slice(0, 18) # slice off the version number from the genes
        adata = anndata.AnnData(X=scipy.io.mmread(mat).tocsr().T, obs=gene, var=bcs)

    del bcs, gene

    print(adata)

    return adata


def basic_process(A):
    '''
        Returns processed adata. obs["counts"] is the column sum
    '''

    adata = A.copy()

    adata.var['counts'] = np.array(adata.X.sum(0))[0]

    adata.var['ngenes'] = np.array((adata.X > 0).sum(0))[0]
    adata = adata[:,adata.var['counts'] > 0]

    adata.layers['log1p'] = np.log1p(adata.X)

    # adata.obs['log10counts'] = np.log10(adata.obs['counts'])

    print(adata)

    return adata

def filter_adata(A, B, obs=False, var=True, by_C=False):
    '''
        Ensures adata_1 and adata_2 have the same var and obs (optional). Filters by the index.
    '''

    adata_1 = A.copy()
    adata_2 = B.copy()

    if by_C:
        adata_1 = adata_1[:, adata_1.var.index.isin(adata_2.var.index)]
        # adata_1 = adata_1[adata_1.obs.index.isin(adata_2.obs.index)]
        return adata_1, adata_2
    if obs:
        adata_1 = adata_1[adata_1.obs.index.isin(adata_2.obs.index)]
        adata_2 = adata_2[adata_2.obs.index.isin(adata_1.obs.index)]
    if var:
        adata_1 = adata_1[:, adata_1.var.index.isin(adata_2.var.index)]
        adata_2 = adata_2[:, adata_2.var.index.isin(adata_1.var.index)]


    return adata_1, adata_2

# Correlations
def _sparse_M_std(X):
    n = X.shape[0]
    return np.sqrt(n * X.multiply(X).sum(0) - np.multiply(X.sum(0), X.sum(0)))

def sparse_M_corr(X, Y):
    '''
        Computes Pearson correlation between X and Y (both in sparse format). Must be same shape.
        X: A_raw[common_obs.index].layers['log1p'] # raw
        Y: B_raw[common_obs.index].layers['log1p']# raw

        X: A.layers['log1p'] # filtered
        Y: B.layers['log1p'] # filtered

        Notes: I changed the axis in sum and shape, need to check if right
    '''
    X_std = _sparse_M_std(X)
    Y_std = _sparse_M_std(Y)
    XY_std = np.multiply(X_std, Y_std)
    n = X.shape[0]
    XY_cov = n*X.multiply(Y).sum(0) - np.multiply(X.sum(0), Y.sum(0))
    R = np.divide(XY_cov, XY_std)
    return np.squeeze(np.asarray(R))

def l1_dist(X, Y):
    '''
        Computes l1 metric between X and Y. May need to modify since I swapped the rows and columns.
        X: A.layers['log1p']
        Y: B.layers['log1p']
    '''
    dist_AA = manhattan_distances(X, X)
    dist_AB = manhattan_distances(X, Y)

    # nkc are the kallisto-cellranger distances
    dist_AB = np.diagonal(dist_AB)

    # ncc are the kallisto-kallisto distances
    AA = []
    for row in dist_AA:
        val = np.partition(row, 1)[1]
        AA.append(val)
    dist_AA = AA

    return dist_AA, dist_AB

def MA(A, B):
    '''
        Computes MA for MA plot
        X: A.var["gene_count"]
        Y: B.var["gene_count"]
    '''
    X = np.array(A.X.mean(axis=1))[:,0]
    Y = np.array(B.X.mean(axis=1))[:,0]

    M_AB = np.log2(X + 1) - np.log2(Y + 1)
    A_AB = 0.5*(np.log2(X + 1) + np.log2(Y + 1))
    return A_AB, M_AB

def barcode_sets(A, B):
    '''
        Computes all, common, and individual barcode sets for two adatas A and B
    '''
    joint = A.var.join(B.var, how = 'outer', lsuffix='-A', rsuffix='-B')
    #joint = joint.fillna(0)

    common = A.var.join(B.var, how = 'inner', lsuffix='-A', rsuffix='-B')

    A_var = A.var.join(B.var, how = 'left', lsuffix='-A', rsuffix='-B')
    A_var = A_var.sort_values(by=['counts-A'], ascending = False)
    A_var = A_var[["counts-A", "ngenes-A"]].rename(columns={"counts-A":"counts", "ngenes-A":"ngenes"})

    B_var = B.var.join(A.var, how = 'left', lsuffix='-B', rsuffix='-A')
    B_var = B_var.sort_values('counts-A', ascending = False)
    B_var = B_var[["counts-B", "ngenes-B"]].rename(columns={"counts-B":"counts", "ngenes-B":"ngenes"})

    return joint, common, A_var, B_var

def compute_tsvd(A):

    adata = A.copy()

    tsvd = TruncatedSVD(n_components=10)
    TSVD = tsvd.fit_transform(adata.layers['log1p'].T)
    # print(TSVD.shape)
    adata.varm['TSVD'] = TSVD

    return adata

def compute_tsne(A):
    adata = A.copy()

    #tsne = TSNE(perplexity=30, metric="euclidean", callbacks=openTSNE.callbacks.ErrorLogger(),n_jobs=8, random_state=42, n_iter=750 )
    tsne = TSNE(perplexity=30, metric="euclidean", n_jobs=10, random_state=42, n_iter=750 )
    adata.varm['TSNE10'] = tsne.fit_transform(adata.varm['TSVD'])

    return adata

def T_sparse_M_std(X):
    n = X.shape[1]
    return np.sqrt(n * X.multiply(X).sum(1) - np.multiply(X.sum(1), X.sum(1)))
def T_sparse_M_corr(X,Y):
    X_std = T_sparse_M_std(X)
    Y_std = T_sparse_M_std(Y)
    XY_std = np.multiply(X_std, Y_std)
    n = X.shape[1]
    XY_cov = n* X.multiply(Y).sum(1) - np.multiply(X.sum(1), Y.sum(1))
    R = np.divide(XY_cov, XY_std)
    return np.squeeze(np.asarray(R))

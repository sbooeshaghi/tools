#!/usr/bin/env python
# coding: utf-8

# # kast: kallisto bus, salmon alevin, star, tenx 

import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import time
import copy
from sklearn.preprocessing import LabelEncoder
from scipy import sparse
import scipy
import anndata
from matplotlib.pyplot import figure
from sklearn.decomposition import TruncatedSVD
import sklearn
from sklearn.metrics import confusion_matrix
import umap
import plotly.plotly as py
import plotly.graph_objs as go
import plotly 
from vpolo.alevin import parser
import anndata
import time
import seaborn as sns
from openTSNE import TSNE
import openTSNE
from openTSNE.callbacks import ErrorLogger
from anndata import AnnData
from sys import argv

script, dataset_shortname, dataset_name = argv

sns.set(style="ticks", context="talk")
plt.style.use("dark_background")



print ("The script is called:", script)
print ("dataset_shortname variable is:", dataset_shortname)
print ("dataset_name variable is:", dataset_name)


# python3.6 kabus_summary_maker.py hgmm1k_v2 '1k 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells (v3 chemistry)' >> log2.txt
# ## Load cellranger and kallisto matrices and gene name to ensembl map
#dataset_name = '1k 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells (v3 chemistry)'
#dataset_shortname = 'hgmm1k_v2'

#dataset_name =  ' 1k Heart Cells from an E18 mouse (v2 chemistry)'
#dataset_shortname = 'heart1k_v2'


# ## Load cellranger and kallisto matrices and gene name to ensembl map


kallisto_folder = '/home/munfred/single_cell_analysis/kallisto_out_single/kallisto_' + dataset_shortname
tenx_folder = '/home/munfred/single_cell_analysis/cellranger_out/cellranger3_' + dataset_shortname +'_out/outs/filtered_feature_bc_matrix'
raw_folder = '/home/munfred/single_cell_analysis/cellranger_out/cellranger3_' + dataset_shortname +'_out/outs/raw_feature_bc_matrix'
star_folder = '/home/munfred/single_cell_analysis/star_out/star_' + dataset_shortname +'_out'
salmon_folder = '/home/munfred/single_cell_analysis/salmon_out/salmon_'+ dataset_shortname +'_out/'


if dataset_shortname == 'pbmc1k_v3': kallisto_folder = '/home/munfred/single_cell_analysis/kallisto_out_single/kallisto_pbmc_1k_v3' 
if dataset_shortname == 'pbmc10k_v3': kallisto_folder = '/home/munfred/single_cell_analysis/kallisto_out_single/kallisto_pbmc_10k_v3' 
if dataset_shortname == 'neuron10k_v3': kallisto_folder = '/home/munfred/single_cell_analysis/kallisto_out_single/kallisto_neuron_10k_v3' 


# ## Load alevin
start = time.time()
# load alevin data...takes a long time
alevin_df = parser.read_quants_bin('/home/munfred/single_cell_analysis/salmon_out/salmon_heart1k_v2_out/')
print('Alevin data shape:', alevin_df.shape)

end = time.time()
totaltime = (end-start)/60
print('Minutes taken to load alevin data:', totaltime)

# load 10x on anndata as sparse csr matrix
tenx_filtered = anndata.AnnData(scipy.io.mmread(os.path.join(tenx_folder,'matrix.mtx.gz')).tocsr().T)
tenxs_barcodes = pd.read_csv(os.path.join(tenx_folder,'barcodes.tsv.gz'), index_col = 0, header = None, names = ['barcode'])
tenxs_barcodes.index = tenxs_barcodes.index.str.slice(0,16,1)
tenx_filtered.obs= tenxs_barcodes
tenx_filtered.var = pd.read_csv(os.path.join(tenx_folder,'features.tsv.gz'), header = None, index_col = 0, names =['ensembl_id', 'gene_name', 'kind'], sep = '\t')
#print('Loaded tenx filtered mtx:\n',tenx_filtered)

alevin_raw = anndata.AnnData(alevin_df)
alevin_raw.obs = pd.DataFrame(index = alevin_df.index.values).rename_axis('barcode')
alevin_raw.var = pd.DataFrame(index = alevin_df.columns.values).rename_axis('ensembl_id')
print('Loaded alevin raw mtx:\n',alevin_raw)
alevin_filtered = alevin_raw[tenx_filtered.obs.index.values]
alevin_filtered.X = scipy.sparse.csr_matrix(alevin_filtered.X)

## load star on anndata
star_raw = anndata.AnnData(scipy.io.mmread(os.path.join(star_folder,'matrix.mtx')).tocsr().T)
star_raw.obs= pd.read_csv(os.path.join(star_folder,'barcodes.tsv'), index_col = 0, header = None, names = ['barcode'])
star_raw.var = pd.read_csv(os.path.join(star_folder,'genes.tsv'), header = None, index_col = 0, names =['ensembl_id','gene_name '], sep = '\t')
print('Loaded star raw mtx:\n',star_raw)
star_filtered = star_raw[tenx_filtered.obs.index.values]

## load kallisto on anndata as sparse crs matrix
kallisto_raw = anndata.AnnData(scipy.io.mmread(os.path.join(kallisto_folder,'genes.mtx')).tocsr())
kallisto_raw.obs= pd.read_csv(os.path.join(kallisto_folder,'genes.barcodes.txt'), index_col = 0, header = None, names = ['barcode'])
kallisto_raw.var = pd.read_csv(os.path.join(kallisto_folder,'genes.genes.txt'), header = None, index_col = 0, names =['ensembl_id'], sep = '\t')
#print('Loaded kallisto raw mtx:\n',kallisto_raw)
kallisto_filtered = kallisto_raw[tenx_filtered.obs.index.values]

#make alldata. its a pun on adata.
alldata =  AnnData.concatenate(tenx_filtered,kallisto_filtered,alevin_filtered,star_filtered, join='outer', batch_categories=['tenx','kallisto','alevin','star'],index_unique='-')
alldata.obs['orig_bc'] = alldata.obs.index.str.split('-').str.get(0)
# check the barcodes in all batches are the same in the same order
if np.array_equal(alldata[alldata.obs.query('batch == "tenx"').index].obs.head()['orig_bc'].values,
               alldata[alldata.obs.query('batch == "kallisto"').index].obs.head()['orig_bc'].values):
    print ('kallisto barcodes are good')
if np.array_equal(alldata[alldata.obs.query('batch == "tenx"').index].obs.head()['orig_bc'].values,
               alldata[alldata.obs.query('batch == "alevin"').index].obs.head()['orig_bc'].values):
    print ('alevin barcodes are good')

if np.array_equal(alldata[alldata.obs.query('batch == "tenx"').index].obs.head()['orig_bc'].values,
               alldata[alldata.obs.query('batch == "star"').index].obs.head()['orig_bc'].values):
    print ('star barcodes are good')

    
#add layer with logged counts for convenience    
alldata.layers['log1p'] = np.log1p(alldata.X)
alldata.obs['counts'] = alldata.X.sum(1)
alldata.obs['log10counts']=np.log10(alldata.obs['counts'])


# make views of each pipeline data for convenience
alevin = alldata[alldata.obs.query('batch == "alevin"').index]
kallisto = alldata[alldata.obs.query('batch == "kallisto"').index]
star = alldata[alldata.obs.query('batch == "star"').index]
tenx = alldata[alldata.obs.query('batch == "tenx"').index]


# # Correlation of cells processed by kallisto and cellranger as a function of expression

corrdf = pd.DataFrame()

for cell in range(len(tenx)):
    ktcorr = np.corrcoef(np.log1p(kallisto.X[cell]).toarray(), np.log1p(tenx.X[cell]).toarray())[0, 1]
    kacorr = np.corrcoef(np.log1p(kallisto.X[cell]).toarray(), np.log1p(alevin.X[cell]).toarray())[0, 1]
    kscorr = np.corrcoef(np.log1p(kallisto.X[cell]).toarray(), np.log1p(star.X[cell]).toarray())[0, 1]
    tacorr = np.corrcoef(np.log1p(tenx.X[cell]).toarray(),     np.log1p(alevin.X[cell]).toarray())[0, 1]
    tscorr = np.corrcoef(np.log1p(tenx.X[cell]).toarray(),     np.log1p(star.X[cell]).toarray())[0, 1]
    sacorr = np.corrcoef(np.log1p(star.X[cell]).toarray(),     np.log1p(alevin.X[cell]).toarray())[0, 1]

    counts = tenx.X[cell].sum()
    corrdf = corrdf.append({
                            'ktcorr':ktcorr,
                            'kacorr':kacorr,
                            'kscorr':kscorr,
                            'tacorr':tacorr,
                            'tscorr':tscorr,
                            'sacorr':sacorr,
                            'counts':counts}, 
                           ignore_index=True)

# # Plot embeddings

# # TSVD

# Do TSVD and keep 50 dimensions for t-SNE and UMAP

tsvd50 = TruncatedSVD(n_components=50)
TSVD50 = tsvd50.fit_transform(alldata.layers['log1p'])
alldata.obsm['TSVD'] = TSVD50
alldata.obsm['TSVD']
print('TSVD variance ratio:\n', tsvd50.explained_variance_ratio_)


# ## UMAP

# Do UAMP ON tsvd 50

alldata.obsm['UMAP'] = umap.UMAP().fit_transform(alldata.obsm['TSVD'])

print('UMAP done.')


# # TSNE-TSVD50

tsne = TSNE(perplexity=30, metric="euclidean", callbacks=openTSNE.callbacks.ErrorLogger(),n_jobs=48, random_state=42 )


# # TSNE-FULL

tsne = TSNE(perplexity=30, metric="euclidean", callbacks=openTSNE.callbacks.ErrorLogger(),n_jobs=24, initialization="random", random_state=42 )

#save the data!!!
alldata.write(os.path.join('./anndatas/' + dataset_shortname + '.h5ad'))

# # Make summary plot composition

def pinplotax(ax, embedding, dim0 = 0, dim1 =1, alldata = alldata, dataset_name=dataset_name, refbatch_name='tenx', dataset_shortname=dataset_shortname):
    refbatch =  alldata[alldata.obs.query(f'batch == "{refbatch_name}"').index]
    refcounts= refbatch.obs['log10counts'].values
    refimbed = refbatch.obsm[embedding]

    colordict= {'kallisto':'lime', 'alevin':'red', 'tenx':'blue','star':'yellow'}

    #fig, ax = plt.subplots(figsize=(20,15))


    for batch in np.unique(alldata.obs['batch']):
        otherimbed = alldata[alldata.obs.query(f'batch == "{batch}"').index].obsm[embedding]

        ax.plot([refimbed[:,dim0][0], otherimbed[:,dim0][0]],
                 [refimbed[:,dim1][0], otherimbed[:,dim1][0]],
                 color = colordict[batch],linewidth=0.4, label = batch )    
        for cellnum in range(1, len(refbatch)):
            ax.plot([refimbed[:,dim0][cellnum], otherimbed[:,dim0][cellnum]],
                     [refimbed[:,dim1][cellnum], otherimbed[:,dim1][cellnum]],
                     color = colordict[batch],linewidth=0.4)
    sc = ax.scatter(refimbed[:,dim0], refimbed[:,dim1],s =10, label="kallisto", c = refcounts, cmap = 'viridis')

    ax.set_ylabel(str(' Dim '+ str(dim1)), fontsize=11)
    ax.set_xlabel(str( 'Dim '+ str(dim0)), fontsize=11)
    #ax.set_title(str(embedding + ' Embedding of samples processed with kallisto bus, Star, Salmon Alevin and Cell Ranger processed \n \n '+ 
    #                 dataset_name + 
    #                 '\n \n A line connects Cell Ranger barcodes (•) to the corresponding barcodes'))
    #cbar = plt.colorbar(sc)
    
    #cbar.set_label('log 10 of molecule counts for the cell', rotation=270, labelpad=30)
    #plt.legend()
    plt.savefig(str('./svgs/' + embedding + '_' + dataset_shortname + '.svg'), dpi=300)
    #plt.savefig(str('./pngs/' + embedding + '_' + dataset_shortname + '.png'), dpi=300)
    #plt.show()
    return sc


fig, ax = plt.subplots(nrows=3, ncols=2,figsize=(22,25))
fsize= 11

####  KNEE PLOT ##########
ax[0,0].plot( np.sort(np.array(tenx.X.sum(1)), axis=None)[::-1],  color ='blue', linewidth=4, label = 'Cell Ranger')
ax[0,0].plot( np.sort(np.array(kallisto_raw.X.sum(1)), axis=None)[::-1], color ='lime', linewidth=1, label = 'kallisto')
ax[0,0].plot( np.sort(np.array(alevin_raw.X.sum(1)), axis=None)[::-1], color ='red',linewidth=1, label = 'Alevin')
ax[0,0].plot( np.sort(np.array(star_raw.X.sum(1)), axis=None)[::-1], color ='yellow',linewidth=1, label = 'Star')

ax[0,0].set_xscale("log", nonposx='clip')
ax[0,0].set_yscale("log", nonposy='clip')
ax[0,0].set_ylabel('Number of reads',fontsize=fsize)
ax[0,0].set_xlabel('barcode',fontsize=fsize)
ax[0,0].set_title('Number of umis per Barcode', fontweight='bold')
ax[0,0].legend()

#### CORRELATION PLOT ####
#ax[0,1].set_ylim(0,1)
corrdf.plot.scatter(x = 'counts', y = 'ktcorr', s =2,ax=ax[0,1], c = 'green',label = 'kallisto-cellranger ')
corrdf.plot.scatter(x = 'counts', y = 'kacorr', s =2,ax=ax[0,1], c = 'steelblue',label = 'kallisto-alevin ')
corrdf.plot.scatter(x = 'counts', y = 'kscorr', s =2,ax=ax[0,1], c = 'lime',label = 'kallisto-star ')
corrdf.plot.scatter(x = 'counts', y = 'tacorr', s =2,ax=ax[0,1], c = 'red',label = 'alevin-cellranger ')
corrdf.plot.scatter(x = 'counts', y = 'tscorr', s =2,ax=ax[0,1], c = 'gold',label = 'star-cellranger ')
corrdf.plot.scatter(x = 'counts', y = 'sacorr', s =2,ax=ax[0,1], c = 'darkorange',label = 'star-alevin')
ax[0,1].set_xscale('log')
ax[0,1].set_title(str('Total counts vs correlation per barcode'), fontweight='bold'  )
ax[0,1].set_xlabel('Cell Ranger counts',fontsize=fsize)
ax[0,1].set_ylabel('Pearson Correlation',fontsize=fsize)
ax[0,1].legend()

#### EMBEDDING PLOTS #####
sc1 = pinplotax(ax = ax[1,0], embedding='TSVD')
ax[1,0].set_title('Truncated SVD first and second components', fontweight='bold')

sc2 = pinplotax(ax = ax[1,1], embedding='UMAP')
ax[1,1].set_title('UMAP on TSVD 50 components', fontweight='bold')

sc3 = pinplotax(ax = ax[2,0], embedding='TSNE-TSVD50')
ax[2,0].set_title('t-SNE on TSVD 50 components',fontweight='bold')

sc4 = pinplotax(ax = ax[2,1], embedding='TSNE-FULL')
ax[2,1].set_title('t-SNE on full count matrices', fontweight='bold')

fig.subplots_adjust(bottom=0.04)
cbar_ax = fig.add_axes([0.15, 0, 0.7, 0.02])
cbar = fig.colorbar(sc1, cax=cbar_ax, orientation="horizontal")
cbar.set_label('log 10 of Cell Ranger counts for each barcode')
#sc = ax.scatter(embedding_dict['K'+ embedding_kind][:,0], embedding_dict['K'+ embedding_kind][:,1], s =20, label="kallisto", c = klogreads, cmap = 'viridis')
#cbar = fig.colorbar()
#cbar.set_label('log 10 of molecule counts for the cell', rotation=270, labelpad=30)

figtitle = (str(' Summary of samples processed with kallisto bus, Star, Salmon Alevin and Cell Ranger \n \n '+ dataset_name.replace("\\n", "\n") + '\n \n A line connects Cell Ranger barcodes (•) to the corresponding barcodes'))

fig.suptitle(figtitle, fontsize=14, fontweight='bold')
plt.subplots_adjust(left=0.125, bottom=0.1, right=0.9, top=0.9 , wspace=0.2, hspace=0.3)

#plt.savefig(str('./svgs/' + 'SUMMARY' + '_' + dataset_shortname + '.svg'), dpi=300, bbox_inches='tight')
plt.savefig(str('./pngs/' + 'SUMMARY' + '_' + dataset_shortname + '.png'), dpi=300, bbox_inches='tight')



print('SAVED THE PLOT!!!!')
print('SAVED THE PLOT!!!!')
print('SAVED THE PLOT!!!!')
print('SAVED THE PLOT!!!!')





import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd

fsize= 15
A_color = '#FF7F0E'
B_color = '#1F77B4'
dotsize = 10
xmax = 1e6
gridalpha = 0.2

def _lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def _plt_color(lst):
    cols=[]
    for l in lst:
        if l>0.25 or l<-0.25:
            cols.append("red")
        elif l<=0.25 and l>=-0.25:
            cols.append('black')
    return cols


def knee_plot(A, ax=None,light=False, **kwargs):
    '''
        Makes knee plot.
        A: adata
        kwargs: [c = _lighten_color(B_color, 0.5), linewidth=2, alpha=1]
    '''
    if light:
        kwargs["c"] = _lighten_color(kwargs["color"], 0.5)
        del kwargs["color"]
        ranked_umi = np.sort(np.array(A.X.sum(0)), axis=None)[::-1]
        return ax.plot(ranked_umi, np.arange(len(ranked_umi)), **kwargs)
    #del kwargs["light"]

    ax = ax or plt.gca()

    if kwargs["label"] == "kallisto":
        ranked_umi = np.sort(np.array(A.X.sum(0)), axis=None)[::-1]
        return ax.plot(ranked_umi[0:A.X.shape[1]], np.arange(A.X.shape[1]), **kwargs)

    ranked_umi = np.sort(np.array(A.X.sum(0)), axis=None)[::-1]
    return ax.plot(ranked_umi, np.arange(len(ranked_umi)), **kwargs)

def knee_plot_settings(A, ax=None):
    ax = ax or plt.gca()

    ranked_umi = np.sort(np.array(A.var["counts"].values), axis=None)[::-1]
    ax.set_xscale('log')
    ax.set_xlim(1,xmax)
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlabel('kallisto UMI counts',fontsize=fsize)
    ax.set_ylabel('Cumulative number of barcodes',fontsize=fsize)
    ax.set_title('',loc='center')
    ax.set_title('A', fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    ax.axhline(y=np.shape(A.X)[1],linewidth=2, color='black', linestyle='--')
    ax.axvline(x=ranked_umi[np.shape(A.X)[1]-1],linewidth=2, color='black', linestyle='--')

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
    return ax

def barcode_ratio_plot(A, B, joint, ax=None, **kwargs):
    ax = ax or plt.gca()

    ax.plot(np.geomspace(1,10e5,100),np.geomspace(1,10e5,100),'gray',linewidth=1 ) # identity line
    ax.scatter(joint['counts-A'].values, joint['counts-B'].values, color ='lightgray', s=dotsize, alpha=0.3, edgecolors = 'none', label = 'Discarded barcodes')
    ax.scatter(A.var['counts'].values, B.var['counts'].values, color ='black', s=dotsize, alpha=0.3, edgecolors = 'none', label = 'Retained barcodes')
    return ax

def barcode_ratio_plot_settings(ax=None):
    ax = ax or plt.gca()

    ax.set_xscale('log')
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlabel('kallisto UMI counts',fontsize=fsize)
    ax.set_ylabel('kallisto velocity UMI counts',fontsize=fsize)
    ax.set_title('',loc='center')
    ax.set_xlim(1,xmax)
    ax.set_ylim(1,xmax)
    ax.set_title('B', fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])

    return ax

def genes_detected_plot(A, ax=None, light=False, **kwargs):
    ax = ax or plt.gca()

    if light:
        kwargs["c"] = _lighten_color(kwargs["c"], 0.3)

    return ax.scatter(A['counts'], A['ngenes'], **kwargs)

def genes_detected_plot_settings(ax=None):
    ax = ax or plt.gca()

    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = 0.5)
    ax.set_xscale('log')
    ax.set_xlim(1,xmax)
    ax.set_yscale("log", nonposy='clip')
    ax.set_ylabel('Genes detected',fontsize=fsize)
    ax.set_xlabel('kallisto UMI counts',fontsize=fsize)
    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    ax.set_title('C', fontweight='bold', fontsize = fsize, loc = 'left' )
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])

    return ax

def cc_raw_plot(common, cc_raw, ax=None, **kwargs):
    ax = ax or plt.gca()

    ax.scatter(x = common["counts-A"], y = cc_raw, **kwargs)
    return ax

def cc_filtered_plot(A, cc_filtered, ax=None, **kwargs):
    ax = ax or plt.gca()

    ax.scatter(x = A.var['counts'], y = cc_filtered, **kwargs)
    return ax

def cc_plot_settings(ax=None):
    ax = ax or plt.gca()

    ax.set_xscale('log')
    ax.set_xlim(1,xmax)
    ax.set_ylim(0,1)
    ax.set_title('D', fontweight='bold', fontsize = fsize, loc = 'left')
    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    ax.set_xlabel('kallisto UMI counts', fontsize = fsize)
    ax.set_ylabel('Pearson Correlation', fontsize = fsize)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])
    return ax

def l1_plot(dist_AA, dist_AB, ax=None, **kwargs):
    ax = ax or plt.gca()

    hist, concat_bins = np.histogram(np.concatenate((dist_AA,dist_AB)), bins='auto')
    hist, ck_bins =  np.histogram(dist_AB, bins='auto')
    hist, cc_bins =  np.histogram(dist_AA, bins='auto')
    best_bins = max([ck_bins,concat_bins,cc_bins], key=len)

    ax.hist(x=dist_AB, bins=best_bins, alpha=0.5, color = B_color, label="Closest Cell Ranger barcode")
    ax.hist(x=dist_AA, bins=best_bins, alpha=0.5, color = A_color, label="Closest kallisto barcode")
    return ax

def l1_plot_settings(ax=None):
    ax = ax or plt.gca()

    ax.set_xlabel('$\ell_1$ Distance', fontsize = fsize)
    ax.set_ylabel('Barcode counts', fontsize = fsize)
    ax.set_title('E', fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1])

    return ax

def tsne_plot(A, ax=None, **kwargs):
    ax = ax or plt.gca()

    ax.scatter(A.varm['TSNE10'][:,0], A.varm['TSNE10'][:,1], **kwargs)

    return ax

def tsne_plot_settings(title, ax=None):
    ax = ax or plt.gca()

    ax.set_ylabel(str(' t-SNE 2'), fontsize=fsize)
    ax.set_xlabel(str( 't-SNE 1'), fontsize=fsize)
    # loc.set_title('t-SNE on TSVD 10 components',fontweight='bold')
    ax.set_title(title, fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.legend()

    return ax

def MA_plot(A_AB, M_AB, ax=None, **kwargs):
    ax = ax or plt.gca()

    cols = _plt_color(M_AB)
    ax.scatter(A_AB, M_AB, alpha=0.05, c=cols)

    return ax

def MA_plot_settings(M_AB, ax=None):
    ax = ax or plt.gca()

    ax.set_ylabel("log2 Fold change", fontsize = fsize)
    ax.set_xlabel("log2 Average gene count", fontsize = fsize)
    ax.set_title('G', fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.set_ylim(-5, 5)
    A_patch = mpatches.Patch(color=A_color, label="kallisto")
    B_patch = mpatches.Patch(color=B_color, label="Cell Ranger")
    same = mpatches.Patch(color='black', label='log2 FC $\leq$ 0.25 ({:.3f})'.format(M_AB[M_AB<=0.25].shape[0]/M_AB.shape[0]))
    ax.arrow(0, 1, 0, 1.5, length_includes_head=True, width=.05, color=A_color)
    ax.arrow(0, -1, 0, -1.5, length_includes_head=True, width=.05, color=B_color)
    ax.legend(handles=[A_patch, B_patch, same])
    return ax

def DE_plot(folder, name, ax=None):
    ax = ax or plt.gca()

    df = pd.read_csv(folder + name + '.csv')
    fold_change_change_colors = {
    '(4,5]':_lighten_color(A_color, 1.4),
    '(3,4]':_lighten_color(A_color, 1.1),
    '(2,3]':_lighten_color(A_color, 0.8),
    '(1,2]':_lighten_color(A_color, 0.5),
    '(0,1]':_lighten_color(A_color, 0.2),
    '(-1,0]':_lighten_color(B_color, 0.2),
    '(-2,-1]':_lighten_color(B_color, 0.5),
    '(-3,-2]':_lighten_color(B_color, 0.8),
    '(-4,-3]':_lighten_color(B_color, 1.1),
    '(-5,-4]':_lighten_color(B_color, 1.4),
    }
    # For some datasets there are no DE genes, so we need this check to just write a text and make not plot
    if len(df)==0:
        ax.text(0.5*(1), 0.5*(1), 'No significant gene sets found ',
            horizontalalignment='center',
            verticalalignment='center',
            fontsize=20, color='black',
            transform=ax.transAxes)
        ax.axis('off')

    # If there are DE genes, then we make a plot
    if len(df)>0:
        ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = 0.5)

        genesets = np.sort(df.gene_set.unique())
        changes = df.change.unique()
        changes = list(reversed(changes))
        changes_dict = {}
        for change_interval in changes:
            changes_dict[change_interval] = []
            for geneset in genesets:
                ngenes_found = np.sum(df.loc[df['gene_set']==geneset]['change']==change_interval)
                changes_dict[change_interval].append(ngenes_found)
        ind = [x for x, _ in enumerate(genesets)]
        bottom_list = [0]*len(genesets)
        for change_interval in changes_dict:
            plt.bar(ind, changes_dict[change_interval], bottom = bottom_list, label = change_interval,
                    alpha = 0.5, width=0.8, color = fold_change_change_colors[change_interval])
            bottom_list = np.array(bottom_list) + np.array(changes_dict[change_interval])

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title='Fold change')

        ax.set_xticks(ticks=ind)
        ax.set_xticklabels( labels=genesets)
        ax.set_ylabel("Number of genes", fontsize=fsize)
        ax.set_xlabel("Gene set", fontsize=fsize)
        for label in ax.get_xmajorticklabels():
            label.set_rotation(30)
            label.set_horizontalalignment("right")

    ax.set_title('H', fontweight='bold', fontsize = fsize, loc = 'left' )

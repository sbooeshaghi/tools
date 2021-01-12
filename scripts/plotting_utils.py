import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable



fsize= 20
A_color = '#FF7F0E'
B_color = '#1F77B4'
dotsize = 10
xmax = 1e6
gridalpha = 0.2
mscale = 4.0

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
    ax.legend(handles[::-1], labels[::-1], fontsize=fsize, markerscale=mscale, loc="upper right")
    return ax

def barcode_ratio_plot(A, B, joint, ax=None, **kwargs):
    ax = ax or plt.gca()

    ax.plot(np.geomspace(1,10e5,100),np.geomspace(1,10e5,100),'gray',linewidth=1 ) # identity line
    ax.scatter(joint['counts-A'].values, joint['counts-B'].values, color ='darkgray', s=dotsize, alpha=0.3, edgecolors = 'none', label = 'Discarded barcodes')
    ax.scatter(A.var['counts'].values, B.var['counts'].values, color ='black', s=dotsize, alpha=0.3, edgecolors = 'none', label = 'Retained barcodes')
    return ax

def barcode_ratio_plot_settings(ax=None):
    ax = ax or plt.gca()

    ax.set_xscale('log')
    ax.set_yscale("log", nonposy='clip')
    ax.set_xlabel('kallisto UMI counts',fontsize=fsize)
    ax.set_ylabel('Cell Ranger UMI counts',fontsize=fsize)
    ax.set_title('',loc='center')
    ax.set_xlim(1,xmax)
    ax.set_ylim(1,xmax)
    ax.set_title('B', fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = gridalpha)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], fontsize=fsize, markerscale=mscale, loc="upper left")

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
    ax.legend(handles[::-1], labels[::-1], fontsize=fsize, markerscale=mscale, loc="lower right")

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
    ax.legend(handles[::-1], labels[::-1], fontsize=fsize, markerscale=mscale, loc="lower right")
    return ax

def l1_plot(dist_AA, dist_AB, ax=None, **kwargs):
    ax = ax or plt.gca()

    hist, concat_bins = np.histogram(np.concatenate((dist_AA,dist_AB)), bins='auto')
    hist, ck_bins =  np.histogram(dist_AB, bins='auto')
    hist, cc_bins =  np.histogram(dist_AA, bins='auto')
    best_bins = min([ck_bins,concat_bins,cc_bins], key=len) # may need to change to max

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
    ax.legend(handles[::-1], labels[::-1], loc="upper right")

    return ax

def mod_l1_plot(counts, nkk, nck, ax_left, ax_right):
    nullfmt = NullFormatter()
    ax_right.yaxis.set_major_formatter(nullfmt)
    ax_left.scatter(counts, nkk, color=A_color, alpha=gridalpha)
    ax_left.scatter(counts, nck, color=B_color, alpha=gridalpha)
    binwidth = 0.25
    xymax = np.max([np.max(np.fabs(counts)), np.max(np.fabs(nkk))])
    lim = (int(xymax/binwidth) + 1) * binwidth

    bins = np.arange(0, lim + binwidth, binwidth)

    hist, concat_bins = np.histogram(np.concatenate((nkk,nck)), bins='auto')
    hist, ck_bins =  np.histogram(nck, bins='auto')
    hist, cc_bins =  np.histogram(nkk, bins='auto')
    best_bins = min([ck_bins,concat_bins,cc_bins], key=len) # may need to change to max




    ax_right.hist(nkk, bins=best_bins, orientation='horizontal', color=A_color, label="kallisto", alpha=0.5)
    ax_right.hist(nck, bins=best_bins, orientation='horizontal', color=B_color, label="Cell Ranger", alpha=0.5)
    ax_right.legend(fontsize=fsize, loc="upper right")

    ax_right.set_ylim(ax_left.get_ylim())
    ax_right.set_xlabel("Barcode counts", fontsize=fsize)
    ax_left.set_xlabel("kallisto UMI counts", fontsize=fsize)
    ax_left.set_ylabel("$\ell_1$ Distance", fontsize=fsize)

    ax_left.set_title("E.1", fontweight='bold', fontsize = fsize, loc = 'left' )
    ax_right.set_title("E.2", fontweight='bold', fontsize = fsize, loc = 'left' )

    return ax_left, ax_right
def counts_l1_dist(counts, nkk, nck, ax):
    nullfmt = NullFormatter()
    ax.scatter(counts, nkk, color=A_color, alpha=gridalpha, label="kallisto", s=dotsize)
    ax.scatter(counts, nck, color=B_color, alpha=gridalpha, label="Cell Ranger", s=dotsize)
    ax.set_xlabel("kallisto UMI counts", fontsize=fsize)
    ax.set_ylabel("$\ell_1$ Distance", fontsize=fsize)

    ax.set_title("E.2", fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.legend(fontsize=fsize, loc="upper left", markerscale=mscale)

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
    ax.legend(fontsize=fsize, markerscale=mscale, loc="upper left")

    return ax

def MA_plot(A_AB, M_AB, ax=None, **kwargs):
    ax = ax or plt.gca()

    cols = _plt_color(M_AB)
    ax.scatter(A_AB, M_AB, alpha=0.1, c=cols)

    return ax

def MA_plot_settings(M_AB, ax=None):
    ax = ax or plt.gca()

    ax.set_ylabel("log$_2$(Fold change)", fontsize = fsize)
    ax.set_xlabel("log$_2$(Average gene count)", fontsize = fsize)
    ax.set_title('G.1', fontweight='bold', fontsize = fsize, loc = 'left' )
    ax.set_ylim(-5, 5)
    A_patch = mpatches.Patch(color=A_color, label="kallisto")
    B_patch = mpatches.Patch(color=B_color, label="Cell Ranger")
    same = mpatches.Patch(color='white', label='log$_2$(FC)$\leq$ 0.25 ({:.3f})'.format(M_AB[M_AB<=0.25].shape[0]/M_AB.shape[0]))
    ax.arrow(0, 1, 0, 1.5, length_includes_head=True, width=.05, color=A_color)
    ax.arrow(0, -1, 0, -1.5, length_includes_head=True, width=.05, color=B_color)
    ax.legend(handles=[A_patch, B_patch, same], fontsize=fsize, loc="upper right")
    return ax

# def DE_plot(folder, name, ax=None):
#     ax = ax or plt.gca()
#
#     df = pd.read_csv(folder + name + '.csv')
#     fold_change_change_colors = {
#     '(4,5]':_lighten_color(A_color, 1.4),
#     '(3,4]':_lighten_color(A_color, 1.1),
#     '(2,3]':_lighten_color(A_color, 0.8),
#     '(1,2]':_lighten_color(A_color, 0.5),
#     '(0,1]':_lighten_color(A_color, 0.2),
#     '(-1,0]':_lighten_color(B_color, 0.2),
#     '(-2,-1]':_lighten_color(B_color, 0.5),
#     '(-3,-2]':_lighten_color(B_color, 0.8),
#     '(-4,-3]':_lighten_color(B_color, 1.1),
#     '(-5,-4]':_lighten_color(B_color, 1.4),
#     }
#     # For some datasets there are no DE genes, so we need this check to just write a text and make not plot
#     if len(df)==0:
#         ax.text(0.5*(1), 0.5*(1), 'No significant gene sets found ',
#             horizontalalignment='center',
#             verticalalignment='center',
#             fontsize=20, color='black',
#             transform=ax.transAxes)
#         ax.axis('off')
#
#     # If there are DE genes, then we make a plot
#     if len(df)>0:
#         ax.grid(color='dimgrey', linestyle='-', linewidth=0.5, which="both", alpha = 0.5)
#
#         genesets = np.sort(df.gene_set.unique())
#         changes = df.change.unique()
#         changes = list(reversed(changes))
#         changes_dict = {}
#         for change_interval in changes:
#             changes_dict[change_interval] = []
#             for geneset in genesets:
#                 ngenes_found = np.sum(df.loc[df['gene_set']==geneset]['change']==change_interval)
#                 changes_dict[change_interval].append(ngenes_found)
#         ind = [x for x, _ in enumerate(genesets)]
#         bottom_list = [0]*len(genesets)
#         for change_interval in changes_dict:
#             ax.bar(ind, changes_dict[change_interval], bottom = bottom_list, label = change_interval,
#                     alpha = 0.5, width=0.8, color = fold_change_change_colors[change_interval])
#             bottom_list = np.array(bottom_list) + np.array(changes_dict[change_interval])
#
#         handles, labels = ax.get_legend_handles_labels()
#         ax.legend(handles[::-1], labels[::-1], title='Fold change', fontsize=fsize, title_fontsize=fsize)
#
#         ax.set_xticks(ticks=ind)
#         ax.set_xticklabels( labels=genesets)
#         ax.set_ylabel("Number of genes", fontsize=fsize)
#         ax.set_xlabel("Gene set", fontsize=fsize)
#         for label in ax.get_xmajorticklabels():
#             label.set_rotation(30)
#             label.set_horizontalalignment("right")
#
#     ax.set_title('H', fontweight='bold', fontsize = fsize, loc = 'left' )

def DE_plot(folder, name, ax=None):
    ax = ax or plt.gca()

    df = pd.read_csv(folder + name + '.csv')
    df["gene_set"] = df.gene_set.astype("str")

    try:
        GO_terms = pd.read_csv("/home/sina/projects/bus/validate/reproduce/GO_terms/GO_" + name + ".txt", sep="\t", header=None)
    except:
        GO_terms = []



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

        df = df[~df.gene_set.str.contains('nan')]
        genesets = np.sort(df.gene_set.unique())

        new = []
        for i in genesets:
            new.append(GO_terms.set_index(1).to_dict()[0][i])

        #print(genesets)
        changes = df.change.unique()
        changes = list(reversed(changes))
        changes_dict = {}
        for change_interval in changes:
            changes_dict[change_interval] = []
            for geneset in genesets:
                ngenes_found = np.sum(df.loc[df['gene_set']==geneset]['change']==change_interval)
                changes_dict[change_interval].append(ngenes_found)
        ind = [x for x, _ in enumerate(genesets)]
        ind = new
        bottom_list = [0]*len(genesets)
        for change_interval in changes_dict:
            ax.bar(ind, changes_dict[change_interval], bottom = bottom_list, label = change_interval,
                    alpha = 0.5, width=0.8, color = fold_change_change_colors[change_interval])
            bottom_list = np.array(bottom_list) + np.array(changes_dict[change_interval])

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles[::-1], labels[::-1], title='Fold change', fontsize=fsize, title_fontsize=fsize)


        ax.set_xticks(ticks=new)
        ax.set_xticklabels( labels=new)
        ax.set_ylabel("Number of genes", fontsize=fsize)
        ax.set_xlabel("Gene set", fontsize=fsize)
        for label in ax.get_xmajorticklabels():
            #label.set_rotation(30)
            label.set_horizontalalignment("right")

    ax.set_title('H', fontweight='bold', fontsize = fsize, loc = 'left' )
    return ax



def QQ_hgmm(dataset_shortname, ax):
    df = pd.read_csv("/home/single_cell_analysis/brain_storm/output/gsea_qq/"+ dataset_shortname + ".csv")
    df.ontology = df.ontology.astype("category")

    from sklearn import preprocessing
    le = preprocessing.LabelEncoder()

    from matplotlib.lines import Line2D
    hg = df[df.mapping.str.contains("Hs")]
    mm = df[df.mapping.str.contains("Mm")]
    legend_elements = [Line2D([0], [0], marker='o', color="w",alpha=0.4, label='Human',markerfacecolor='k', markersize=10),
                      Line2D([0], [0], marker='s', color='w', alpha=0.2, label='Mouse', markerfacecolor='grey', markersize=10)]

    #c = le.fit_transform(df.ontology.values)

    c = le.fit_transform(hg.ontology.values)
    scatter = ax.scatter(hg.uniform_log, hg.p_log, c=c)

    c = le.fit_transform(mm.ontology.values)
    ax.scatter(mm.uniform_log, mm.p_log, c=c, marker='s')

    l1 = ax.legend(handles=legend_elements, loc="lower right", title="Species", fontsize=fsize-5, title_fontsize=fsize-5)

    ax.plot(hg.uniform_log, hg.upper, color='k', alpha=0)
    ax.plot(hg.uniform_log, hg.lower, color='k', alpha=0)

    ax.plot(mm.uniform_log, mm.upper, color='grey', alpha=0)
    ax.plot(mm.uniform_log, mm.lower, color='grey', alpha=0)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.fill_between(hg.uniform_log, hg.lower, hg.upper, color='black', alpha='0.4')
    ax.fill_between(mm.uniform_log, mm.lower, mm.upper, color='grey', alpha='0.2')



    # Produce a legend for the ranking (colors). Even though there are 40 different
    # rankings, we only want to show 5 of them in the legend.

    l2 = ax.legend(*(scatter.legend_elements()[0], le.classes_), loc="upper left", title="Ontology", fontsize=fsize-5, title_fontsize=fsize-5)
    ax.add_artist(l1)


    ax.set_xlabel("Expected -log$_{10}$(p)", fontsize=fsize)
    ax.set_ylabel("Observed -log$_{10}$(p)", fontsize=fsize)
    df = df[df.label.astype(str).values != 'nan']
    with open("GO_"+dataset_shortname + ".txt", 'w') as f:


        for i, txt in enumerate(df["rank"]):
            f.write(str(i+1) + "\t" + df.label.values[i] + "\n")
            print(df.label.values[i])
            ax.annotate(txt, (df.uniform_log[i], df.p_log[i]), fontsize=15)
    ax.set_title("G.2", fontweight='bold', fontsize = fsize, loc = 'left' )
    return ax

def QQ_plot(dataset_shortname, ax):
    if "mm" in dataset_shortname:
        return QQ_hgmm(dataset_shortname, ax)

    df = pd.read_csv("/home/single_cell_analysis/brain_storm/output/gsea_qq/"+ dataset_shortname + ".csv")
    df.ontology = df.ontology.astype("category")

    from sklearn import preprocessing
    le = preprocessing.LabelEncoder()


    c = le.fit_transform(df.ontology.values)

    scatter = ax.scatter(df.uniform_log, df.p_log, c=c)


    ax.plot(df.uniform_log, df.upper, color='k', alpha=0)
    ax.plot(df.uniform_log, df.lower, color='k', alpha=0)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.fill_between(df.uniform_log, df.lower, df.upper, color='black', alpha='0.4')



    # Produce a legend for the ranking (colors). Even though there are 40 different
    # rankings, we only want to show 5 of them in the legend.

    l2 = ax.legend(*(scatter.legend_elements()[0], le.classes_), loc="best", title="Ontology", fontsize=fsize-5, title_fontsize=fsize-5)


    ax.set_xlabel("Expected -log$_{10}$(p)", fontsize=fsize)
    ax.set_ylabel("Observed -log$_{10}$(p)", fontsize=fsize)
    df = df[df.label.astype(str).values != 'nan']
    with open("GO_"+dataset_shortname + ".txt", 'w') as f:


        for i, txt in enumerate(df["rank"]):
            f.write(str(i+1) + "\t" + df.label.values[i] + "\n")
            print(df.label.values[i])
            ax.annotate(txt, (df.uniform_log[i], df.p_log[i]), fontsize=15)
    ax.set_title("G.2", fontweight='bold', fontsize = fsize, loc = 'left' )
    return ax

# def QQ_plot(folder, name, ax=None):
#     ax = ax or plt.gca()
#     df = pd.read_csv("/home/single_cell_analysis/brain_storm/output/gsea_qq/"+ name + ".csv")
#     df.ontology = df.ontology.astype("category")
#
#     from sklearn import preprocessing
#     le = preprocessing.LabelEncoder()
#
#
#
#     c = le.fit_transform(df.ontology.values)
#
#     scatter = ax.scatter(df.uniform_log, df.p_log, c=c)
#     ax.plot(df.uniform_log, df.upper, color='k', alpha=0)
#     ax.plot(df.uniform_log, df.lower, color='k', alpha=0)
#
#     lims = [
#         np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
#         np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
#     ]
#
#     # now plot both limits against eachother
#     ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
#     ax.set_aspect('equal')
#     ax.set_xlim(lims)
#     ax.set_ylim(lims)
#
#     ax.fill_between(df.uniform_log, df.lower, df.upper, color='grey', alpha='0.2')
#
#     legend1 = ax.legend(*(scatter.legend_elements()[0], le.classes_), loc="upper left", title="Ontology", fontsize=fsize, title_fontsize=fsize, markerscale=mscale-2)
#
#     ax.set_xlabel("Expected -log$_{10}$(p)", fontsize=fsize)
#     ax.set_ylabel("Observed -log$_{10}$(p)", fontsize=fsize)
#     ax.set_title("G.2", fontweight='bold', fontsize = fsize, loc = 'left' )
#
#     return ax



def make_hist(A, B, orientation="vertical", ax=None):
    hist, concat_bins = np.histogram(np.concatenate((A,B)), bins='auto')
    hist, A_bins =  np.histogram(A, bins='auto')
    hist, B_bins =  np.histogram(B, bins='auto')

    best_bins = min([A_bins,concat_bins,B_bins], key=len) # may need to change to max

    ax.hist(A, bins=best_bins, orientation=orientation, color=A_color, label="kallisto", alpha=1)
    ax.hist(B, bins=best_bins, orientation=orientation, color=B_color, label="Cell Ranger", alpha=1)
    return ax

def make_scatter_hist(dist_AA, dist_BB, dist_AB, dist_BA, ax=None):
    x = dist_AA
    y = dist_AB

    xx = dist_BA
    yy = dist_BB

    # the scatter plot:
    ax.scatter(x, y, label="kallisto", color=A_color)
    ax.scatter(xx, yy, label="Cell Ranger", color=B_color)
    ax.set_aspect(1.)

    # create new axes on the right and on the top of the current axes
    # The first argument of the new_vertical(new_horizontal) method is
    # the height (width) of the axes to be created in inches.
    divider = make_axes_locatable(ax)
    axHistx = divider.append_axes("top", 1.5, pad=0.075, sharex=ax)
    axHisty = divider.append_axes("right", 1.5, pad=0.075, sharey=ax)

    # make some labels invisible
    axHistx.xaxis.set_tick_params(labelbottom=False)
    axHisty.yaxis.set_tick_params(labelleft=False)


    ## Right histogram,  cellranger-cellranger, kallisto-cellranger,
    axHisty = make_hist(y, yy, orientation="horizontal", ax=axHisty)

    # kallisto-kallisto, cellranger-kallisto
    axHistx = make_hist(x, xx, orientation="vertical", ax=axHistx)


    # the xaxis of axHistx and yaxis of axHisty are shared with axScatter,
    # thus there is no need to manually adjust the xlim and ylim of these
    # axis.
    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)

    ax.set_xlabel("$\ell_1$ to nearest kallisto", fontsize=fsize)
    ax.set_ylabel("$\ell_1$ to nearest Cell Ranger", fontsize=fsize)
    axHistx.set_ylabel("Barcode counts", fontsize=fsize-8)
    axHisty.set_xlabel("Barcode counts", fontsize=fsize-8)

    axHistx.set_title("E.1", fontweight='bold', fontsize = fsize, loc = 'left' )

    axHistx.legend(fontsize=fsize-5, loc="upper right")
    return axHistx

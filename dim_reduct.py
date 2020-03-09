#!/usr/bin/python

""" module to run various dimensionality reduction, primarily tSNE and others.
also includes classifier work and support plotting features
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gct, os, math, plottools, gt, sys, mpld3

from sklearn import (manifold, decomposition, ensemble,
                     discriminant_analysis)
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import normalize, scale
from sklearn import svm
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import (MaxAbsScaler, RobustScaler, StandardScaler)

# need checks for dataframe name and samples as rows
def tsne2(df, h, shape='type', outpath='dflt', data=False, legend='full', rm=False, annot=False,
          px='dflt', lr='dflt', test=False, title='dflt', inter=False, cmap=None, labels='dflt',
          hue='name', scaling=None):
    """  !!!!  make sure samples are rows!!!!
    revamped tsne to work with seaborn plotting, pass datafram and header in to generate
    the shape argument can be used to make a header element be a different shape in the plots
     rm is remove features, default is False but can be set to true or set to a specific number
     legend can be 'False', 'brief' or 'full'
     px and lr as 'dflt' will perform a range of both, specific numbers can be entered intstead

     currently relplot is only small, can't figure out how to get it big. but ax works fine large

     shape has max of 8 unique shapes that can be made

    scaling default is set to None, can be 'maxaxs' MaxAbs or 'robust' Robust, 'std' StandardScalar"""

    x, y = ('tsne_x', 'tsne_y')

    if outpath == 'dflt':
        outpath = gt.dflt_outpath(fldr_name='dim_reduct')

    # check if df is transposed or not, make sure samples are rows
    # if ':' in df.columns.values[0]:
    #     df = df.T

    # check for defaults  and assign title
    if isinstance(title, str) and title != 'dflt':
        mytitle = title
    else:
        mytitle = df.name


    if scaling is not None:
        if scaling == 'maxabs':
            df = MaxAbsScaler().fit_transform(df)
        elif scaling == 'robust':
            df = RobustScaler().fit_transform(df)
        elif scaling == 'std':
            df = StandardScaler(with_mean=False).fit_transform(df)
        else:
            print('error with scaling')

    # check if remove features is true, auto removes to 75 by default
    if rm is not False:
        if rm is True:
            # set number of features
            n = 75
        else:
            n = rm
        pca = PCA(n_components=n)
        print('features removed')
        df = pca.fit_transform(df)

    # read in perplexity and lr or assign default ranges
    if isinstance(px, str):
        if px == 'dflt':
            perplexity = [5, 10, 15, 30]
    else:
        if isinstance(px, int):
            perplexity = [int(px)]
        else:
            perplexity = [int(x) for x in px]
    if isinstance(lr, str):
        if lr == 'dflt':
            learning_rate = [10, 100, 250, 400, 600, 800]
    else:
        if isinstance(lr, int):
            learning_rate = [lr]
        else:
            learning_rate = [int(x) for x in lr]

    hdrs = ['n', 'rmf', 'px', 'lr']

    # the default labels assigned are 'label' or 'name'
    if isinstance(labels, str) and labels == 'dflt':
        try:
            labels = h['label']
        except:
            labels = h['name']

    # upstream args
    for p in perplexity:
        for l in learning_rate:
            vals = [df.shape[0], rm, p, l]
            mytitle = ''
            for hdr, v in zip(hdrs, vals):
                mytitle += ('_{}={}'.format(hdr, v))
            tsne = manifold.TSNE(perplexity=p, init='pca', learning_rate=l)
            X = tsne.fit_transform(df)
            h['tsne_x'] = X[:,0]
            h['tsne_y'] = X[:,1]
            # for hdr, v in zip(hdrs, vals):
            #     mytitle += ('_{}={}'.format(hdr, v))
            mypath = os.path.join(outpath, title + '.png')
            if inter is False:
                seaborn_scatter(h, mytitle, outpath, hue=hue, shape=shape, annot=annot, legend=legend, cmap=cmap, ptype='rel')
            elif inter is True:
                fig, ax = seaborn_scatter(h, mytitle, outpath, hue=hue, shape=shape, ptype='ax', save=False, cmap=cmap)
                scatter = ax.scatter(h[x].tolist(), h[y].tolist(), alpha=0.001)
                tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
                mpld3.plugins.connect(fig, tooltip)
                myoutpath = os.path.join(outpath, mytitle + '.html')
                mpld3.save_html(fig, myoutpath)

            if data is True:
                h.to_excel(mypath.replace('.png', '.xlsx'))
            if test is True:
                sys.exit('test mode, quitting after 1')
            plt.close()


# need to test out
def html_scatter(fig, ax, h, x_col, y_col, title, labels='dflt', outpath='dflt'):
    outpath = gt.check_dfltarg(outpath, gt.dflt_outpath('foo'))
    labels = gt.check_dfltarg(labels, h.label)
    scatter = ax.scatter(h[x_col].tolist(), h[y_col].tolist(), alpha=0.001)
    tooltip = mpld3.plugins.PointLabelTooltip(scatter, labels=labels)
    mpld3.plugins.connect(fig, tooltip)
    myoutpath = os.path.join(outpath, title + '.html')
    mpld3.save_html(fig, myoutpath)


def seaborn_scatter(h, title, outpath, shape='type', hue='name', annot=False, legend='brief',
                    save=True, cmap=None, labels='dflt', ptype='ax', x='tsne_x', y='tsne_y',
                    html=False):
    """ uses output of tSNE generation (header file with x and y distance coordinates)
     to generate plot using seaborn instead of get_to_scatter and manual plotting
     legend can be 'full' or 'none' or 'brief'
    can be 'ax' type scatter or alternatively 'rel' for a figure level one
    legends don't work well with ax, but do for relplots. ax required for html """

    if ptype == 'ax':
        fig, ax = plt.subplots()
        fig.set_size_inches(8.5, 8.5)

    # limit a huge legend by setting to 0 if more than 20 long
    try:
        if len(h[hue].unique()) > 25:
            legend = None
    except:
        legend = None

    if isinstance(labels, str) and labels == 'dflt':
        try:
            labels = h['label']
        except:
            labels = h.index.values

    # assume name as color
    kwargs = {'x':x, 'y':y, 'alpha':.8, 'data': h, 'legend':legend}
    # default is 'name'
    kwargs['hue'] = hue

    if cmap is not None:
        kwargs['palette'] = cmap
        sns.set_palette(cmap)

    if shape == 'b' or shape == 'batch':
        kwargs['style'] = 'batch'
    elif shape == 'c' or shape =='cell':
        kwargs['style'] = 'cell'
    elif shape == 't' or shape =='type':
        kwargs['style'] = 'type'
    elif shape == 'p' or shape =='plate':
        kwargs['style'] = 'plate'
    else:
        kwargs['style'] = shape

    try:
        dosecol = h.dose
        numdoses = len(gt.test_only(h).dose.dropna().unique())
    except:
        numdoses = 0
    if 'dose' in h.columns and numdoses == 2:
        kwargs['size'] = 'dose'
        kwargs['sizes'] = (75, 175)
    if 'dose' in h.columns and numdoses >= 3:
        kwargs['size'] = 'dose'
        kwargs['sizes'] = (50, 350)
    else:
        kwargs['s'] = 50

    # plot one of two ways, individual ax object or full g/figure object
    if ptype == 'ax':
        kwargs['ax'] = ax
        kwargs['legend'] = None
        ax = sns.scatterplot(**kwargs)
        ax.set_title(title, fontsize=24)
    else:
        g = sns.relplot(**kwargs)
        g.fig.suptitle(title)

    if annot is True:
        for l, x, y in zip(labels, h[x].values, h[y].values):
            plt.annotate(l, xy=(x, y), xytext=(-2, 2),
                       textcoords='offset points', ha='right', va='bottom')

    mypath = os.path.join(outpath, title + '.png')

    if save is True:
        plt.savefig(mypath)
        plt.close()
    else:
        if ptype == 'ax':
            return fig, ax
        else:
            return g


def clustered_order(df):
    """ just to get the ordering of the samples out of a wardrs linkage association agglomerative
     clustering of data """
    linked = linkage(df.T, 'ward')
    clust_list = leaves_list(linked)
    return clust_list


def make_dendrogram(df, labels='dflt', orient='top', outpath=True, trunc=False, res=False):
    """ uses scipy ward clustering to form dendrogram, trunc folds up last detailed splits """

    linked = linkage(df.T, 'ward')

    fig, ax = plt.subplots()

    if isinstance(labels, str) and labels == 'dflt':
        labels = df.columns

    if trunc is False:
        dend = dendrogram(linked,
                   orientation=orient,
                   labels=labels,
                   distance_sort='descending',
                   show_leaf_counts=True,
                   leaf_rotation=90)

    else:
        dend = dendrogram(linked,
                          orientation=orient,
                          labels=labels,
                          distance_sort='descending',
                          show_leaf_counts=True,
                          truncate_mode='lastp',
                          show_contracted=True,
                          leaf_rotation=90)

    name = gt.dflt_name(df)

    fig.suptitle(f'{name} n={len(df.columns)} dendrogram', fontweight='bold')

    plt.tight_layout()

    if outpath is True:
        outpath = gt.dflt_outpath(fn=name + '_dend.png')
        plt.savefig(outpath)
        plt.close()
    elif outpath != False:
        plt.savefig(outpath)
        plt.close()
    else:
        if res is False:
            return ax
        else:
            return dend


def hier_clustering(df, nclusters):
    """ uses sklearn to agglomerate clusters to group into x clusters, returns list of category identities """

    cluster = AgglomerativeClustering(n_clusters=nclusters, affinity='euclidean', linkage='ward')

    results = cluster.fit_predict(df.T)

    return results


def kmeans_cluster(df, clusters, nc2=None):
    """ conducts KMeans clustering with given dataframe, with provided number of clusters
    and optionally the nc2 (number of clusters 2) which defines the end point of a sweep of
    number of clusters, and the code will then loop through populating a dictionary with cluster
    labels for each column across the range of cluster size
    default behavior is that each ROW is an observation (sample) and columns are features, so the
    dataset is transposed for the fit method """

    if nc2 is not None:
        nclusters = range(clusters, nc2+1)
    else:
        nclusters = [clusters]
    res_dict = {}
    for n in nclusters:
        model = KMeans(n_clusters=n)
        model.fit(df.T)
        res = model.predict(df.T)
        print(f'label length: {len(res)}')
        res_dict[n] = res

    if len(nclusters) == 1:
        return res
    else:
        return res_dict


def plot_kmeans_clusters(df, cat_dict):
    """ plot correlation matrices for a dataset according to the passed label dictionary
    this is designed to be used with the dictionary output of 'kmeans_clusters' above"""

    outpath = gt.dflt_outpath(fn=df.name)

    for cat, labels in cat_dict.items():

        myoutpath = outpath + f'_{cat}_clusters.png'

        plottools.plot_correlation_matrix(df, labels=labels, sort=True, outpath=myoutpath, sparselabel=True)


def svm_classifier(data, classes):
    """ creates a SVM multi-class classifier with passed dataframe and list of classes. returns the
     classifier object """
    X = np.array(data)
    y = classes
    if len(data.index) != len(y):
        X = data.T
    clf = svm.SVC()
    try:
        clf.fit(X, y)
        return clf
    except:
        print(f'problem with data dims {len(data.index)}, {len(data.columns)} and label length {len(y)}')
        return None


def reduce_features(d, met, n):
    """ reduce the feature dimennsion (rows) of passed dataframe. met(hod) can be 'pca' or 'tsvd' with
    then the target number of components"""
    try:
        if met == 'pca':
            pca = PCA(n_components=n)
            dat = pca.fit_transform(d)
        elif met == 'tsvd':
            tsvd = TruncatedSVD(n_components=n)
            dat = tsvd.fit_transform(d)
        else:
            print('no matching method found - "pca" or "tsvd"')
            dat = d
    except ValueError:
        print('!!! error, ran into broadcast error during feat recuction - ', met, n)
        dat = d
    return dat


def run_methods(df, h, outpath='dflt', hdr='dflt', hue='name', mets='dflt', shape=None, labels='dflt', scaling=None):
    """ run a selection of alternate dimension reduction techniques"""
    method_results = dict()

    n_neighbors = 15

    try:
        hdr = gt.check_dfltarg(hdr, df.name)
    except:
        hdr = df.columns[0].split(':')[0]
    outpath = gt.check_dfltarg(outpath, gt.dflt_outpath(fldr_name='dim_reduct'))
    labels = gt.check_dfltarg(labels, h.label)
    mets = gt.check_dfltarg(mets, ['PCA','ISO','MDS','LLE'])

    if ':' in df.columns.values[0]:
        df = df.T

    if scaling is not None:
        if scaling == 'maxabs':
            df = MaxAbsScaler().fit_transform(df)
        elif scaling == 'robust':
            df = RobustScaler().fit_transform(df)
        elif scaling == 'std':
            df = StandardScaler(with_mean=False).fit_transform(df)
        else:
            print('error with scaling')

    if 'PCA' in mets:
        # Projection on to the first 2 principal components
        method_results['PCA'] = decomposition.PCA(n_components=2).fit_transform(df)

    if 'ISO' in mets:
        # Isomap projection of the digits dataset
        method_results['ISO'] = manifold.Isomap(n_neighbors, n_components=2).fit_transform(df)

    if 'MDS' in mets:
        # MDS  embedding of the digits dataset
        method_results['MDS'] = manifold.MDS(n_components=2, n_init=1, max_iter=100).fit_transform(df)

    if 'LLE' in mets:
        # Locally linear embedding
        method_results['LLE'] = manifold.LocallyLinearEmbedding(n_neighbors, n_components=2,
                                                            method='modified').fit_transform(df)

    for met, dat in method_results.items():
        xcol, ycol = f'{met}_x', f'{met}_y'
        h[xcol] = dat[:, 0]
        h[ycol] = dat[:, 1]

        title = hdr + ' ' + met

        seaborn_scatter(h, title, outpath, hue=hue, x=xcol, y=ycol, labels=labels, shape=shape, legend='brief')


def setup_tsne(dat, labels, outpath, hdr):
    """ takes a transposed dataframe (sample rows), labels of colors (should be dict to create
        legend to include with figures) and the folder

        what's hdr?? """

    if ':' in dat.columns.values[0]:
        dat = dat.T

    # tsne permutations
    # norm = [False, True]
    # feat_red_met = ['tsvd', 'pca', 'none']
    # num_feats = [50, 100, 300]

    feat_red_met = ['none']
    num_feats = [50]

    # perplexity = [5, 10, 15, 20, 30]
    # init = ['pca', 'random']
    # learning_rate = [10, 50, 100, 250, 500, 1000]
    perplexity = [5, 10, 15]
    init = ['pca']
    learning_rate = [10, 100, 250]

    hdrs = ['n', 'nor', 'rfm', 'nf', 'px', 'in', 'lr']

    # upstream args
    for rm in feat_red_met:
        for f in num_feats:
            for p in perplexity:
                for i in init:
                    for l in learning_rate:
                        vals = [dat.shape[0], hdr, rm, f, p, i, l]
                        title = 'tsne'
                        for h, v in zip(hdrs, vals):
                            title += ('_{}={}'.format(h, v))
                        run_tsne(dat, rm, f, p, i, l, labels, title, outpath)


def run_tsne(dat, rm, f, p, i, l, labels, title, outpath):
    """ remember rows are samples, columns features => transpose df up front """

    if rm is not 'none':
        dat = reduce_features(dat, rm, f)
    tsne = manifold.TSNE(perplexity=p, init=i, learning_rate=l)

    try:
        X = tsne.fit_transform(dat)
        # make_plot(X, labels, title, outpath)
        df = get_to_scatter(X, labels)
        make_better_plot(df, title, outpath)
        #make_dose_plot(df, title, outpath)
    except AssertionError:
        print('!!!!!  error - ', title)

    # dat = pd.DataFrame(tsne.embedding_, columns=['x','y'])
    # dat.to_csv(outpath + '.txt', sep='\t')


def make_plot(dat, labels, title, outpath):
    """ basic scatter plot """
    print(title)
    fig, ax = plt.subplots()
    ax.scatter(dat[:, 0], dat[:, 1], c=labels)
    ax.set_title(title.replace('_', ' ').replace('=',':'))
    plt.savefig(outpath + '.png')
    plt.close()


def get_to_scatter(da, labels):
    """ the hard-coded method to separate explicitly plotted colors  and labels
    this isn't needed when just using seaborn plots """
    try:
        labels.fillna('blank', inplace=True)
    except:
        pass
    try:
        labels = labels.values
    except:
        pass
    # transform array to dataframe and give coord headers
    df = pd.DataFrame(da, columns=['x', 'y'])
    # assign the full name as the full plate:well:name text
    df['full'] = labels
    # is this necessary?
    # pull the category label from the last ':' delimited part of column name
    try:
        df['category'] = [x.split(':')[-1] for x in labels]
    except:
        df['category'] = labels
    # assign color digits according to category label
    cat_list = sorted(df['category'].unique())
    num_cats = len(cat_list)
    dic = dict()
    # use the presence of decimial point '.' to see if there's doses or not
    # assumed to be in second to last position if so
    if any(['.' in x for x in labels]):
        hdrs = ['category', 'full', 'x', 'y', 'dose', 'color']
        df['dose'] = [x.split(':')[-2] for x in labels]
    else:
        hdrs = ['category', 'full', 'x', 'y', 'color']
    if num_cats < 10:
        cmap = plt.get_cmap('tab10')
        color_coords = np.linspace(0.01, .98, num=10)
        for l, n in zip((cat_list), color_coords[:num_cats]):
            dic[l] = cmap(n)
    if num_cats > 10:
        cmap = plt.get_cmap('nipy_spectral')
        color_coords = np.linspace(0.05, .95, num=num_cats)
        for l, n in zip((cat_list), color_coords):
            dic[l] = cmap(n)
    df['color'] = df['category'].apply(lambda x: dic[x])
    df = df[hdrs]
    return df


def make_better_plot(df, title, outpath):
    """ scatter plot with different colors, assumes df with 'x', 'y', 'category', 'color', 'full' fields """
    print(title)
    subdirs = ['figs', 'figs_annot', 'coords']
    for sd in subdirs:
        try:
            os.mkdir(os.path.join(outpath, sd))
        except:
            pass
    fig, ax = plt.subplots()
    for name, group in df.groupby('category'):
        ax.scatter(group['x'], group['y'], c=group['color'], label=name)
    ax.set_title(title.replace('_', ' ').replace('=',':'))
    #if len(df) < 16:
    #    ax.legend()
    outfold = os.path.join(outpath, 'figs')
    plt.savefig(outpath + 'figs/' + title + '.png')
    for l, x, y in zip(df.index.values, df['x'].values, df['y'].values):
        plt.annotate(l, xy=(x, y), xytext=(-2,2),
                     textcoords='offset points', ha='right', va='bottom')
    outfold = os.path.join(outpath, 'figs_annot')
    plt.savefig(os.path.join(outfold, title + '.png'))
    plt.close()
    df_save = df[['category', 'full', 'x', 'y']]
    outfold = outfold.replace('figs_annot', 'coords')
    df_save.to_excel(os.path.join(outfold, title + '.xlsx'))


def make_dose_plot(df, title, outpath):
    """ scatter plot with different colors, assumes df with 'x', 'y', 'category', 'color', 'full', 'dose' fields """
    print(title)
    doses = sorted(df['dose'].unique())
    sizes = [(doses.index(x) + 1) * 50 for x in doses]
    sizedict = dict(zip(doses, sizes))
    df['size'] = df['dose'].apply(lambda x: sizedict[x])
    fig, ax = plt.subplots()
    for name, group in df.groupby('category'):
        ax.scatter(group['x'], group['y'], c=group['color'], s=group['size'], label=name)
    ax.set_title(title.replace('_', ' ').replace('=',':'))
    #ax.legend()
    plt.savefig(outpath + 'figs/' + title + '.png')
    for l, x, y in zip(df.index.values, df['x'].values, df['y'].values):
        plt.annotate(l, xy=(x, y), xytext=(-2,2),
                     textcoords='offset points', ha='right', va='bottom')
    plt.savefig(outpath + 'figs_annot/' + title + '.png')
    plt.close()
    df_save = df[['category', 'full', 'x', 'y']]
    df_save.to_excel(outpath + 'coords/' + title + '.xlsx')


def round_num(n):
    n = round(float(n), 2)
    return n


def calc_rep_dist(file):
    f = pd.read_excel(file)
    dist = {}
    for n in f.category.unique():
        ps = f[f['category']==n]
        p1 = ps.iloc[0][['x', 'y']].values
        p2 = ps.iloc[1][['x', 'y']].values
        dist[n] = math.sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p1[1])*(p1[1] - p1[1]))
    vals = list(dist.values())
    avg = round_num(np.mean(vals))
    med = round_num(np.median(vals))
    x_range = max(f['x'])
    y_range = max(f['y'])


    return avg, med





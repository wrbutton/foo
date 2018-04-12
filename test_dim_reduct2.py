#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gct, pt, os, math

from sklearn import (manifold, decomposition, ensemble,
                     discriminant_analysis)
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import normalize, scale


def reduce_features(d, met, n):
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


def run_methods(X, outpath, labels, hdr):
    method_results = dict()

    n_neighbors = 15

    # Projection on to the first 2 principal components
    # method_results['X_pca'] = decomposition.TruncatedSVD(n_components=2).fit_transform(X)
    method_results['X_tSVD'] = decomposition.TruncatedSVD(n_components=2).fit_transform(X)
    method_results['X_PCA'] = decomposition.PCA(n_components=2).fit_transform(X)

    # Projection on to the first 2 linear discriminant components
    # method_results['X_lda'] = discriminant_analysis.LinearDiscriminantAnalysis(n_components=2).fit_transform(X, labels)

    # Isomap projection of the digits dataset
    method_results['X_iso'] = manifold.Isomap(n_neighbors, n_components=2).fit_transform(X)

    # MDS  embedding of the digits dataset
    method_results['X_MDS'] = manifold.MDS(n_components=2, n_init=1, max_iter=100).fit_transform(X)

    # Random Trees embedding of the digits dataset
    # X_trans = ensemble.RandomTreesEmbedding(n_estimators=200, random_state=0,
    #                                                          max_depth=5).fit_transform(X)
    # method_results['X_randt'] = decomposition.TruncatedSVD(n_components=2).fit_transform(X_trans)

    # Spectral embedding of the digits dataset
    # method_results['X_se'] = manifold.SpectralEmbedding(n_components=2, random_state=0,
    #                                                   eigen_solver="arpack").fit_transform(X)

    # Locally linear embedding of the digits dataset
    # methods = ['standard', 'ltsa', 'hessian', 'modified']
    methods = ['modified']
    met_names = ['X_lle_std', 'X_lle_ltsa', 'X_lle_hess', 'X_lle_mod']
    for n, m in zip(met_names, methods):
        method_results[n] = manifold.LocallyLinearEmbedding(n_neighbors, n_components=2,
                                                            method=m).fit_transform(X)

    for met, dat in method_results.items():
        title = met.split('_',1)[-1] + '_' + 'n=' + str(X.shape[0]) + '_' + hdr
        # make_plot(dat, labels, title , outpath)
        df = get_to_scatter(dat, labels)
        # make_better_plot(df, title, outpath)
        make_dose_plot(df, title, outpath)



def setup_tsne(dat, labels, path, hdr):
    """ takes a transposed dataframe (sample rows), labels of colors (should be dict to create
        legend to include with figures) and the folder """

    # tsne permutations
    # norm = [False, True]
    # feat_red_met = ['tsvd', 'pca', 'none']
    # num_feats = [50, 100, 300]

    feat_red_met = ['none']
    num_feats = [50]

    # perplexity = [5, 10, 15, 20, 30]
    # init = ['pca', 'random']
    # learning_rate = [10, 50, 100, 250, 500, 1000]
    perplexity = [5, 10, 18, 30]
    init = ['pca']
    learning_rate = [10, 100, 250, 500, 1000]

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
                        run_tsne(dat, rm, f, p, i, l, labels, title, path)


def run_tsne(dat, rm, f, p, i, l, labels, title, outpath):
    """ remember rows are samples, columns features => transpose df up front """

    if rm is not 'none':
        dat = reduce_features(dat, rm, f)
    tsne = manifold.TSNE(perplexity=p, init=i, learning_rate=l)

    try:
        X = tsne.fit_transform(dat)
        # make_plot(X, labels, title, outpath)
        df = get_to_scatter(X, labels)
        # make_better_plot(df, title, outpath)
        make_dose_plot(df, title, outpath)
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



def get_to_scatter(da, full_name_list):
    # transform array to dataframe and give coord headers
    df = pd.DataFrame(da, columns=['x', 'y'])
    # assign the full name as the full plate:well:name text
    df['full'] = full_name_list
    # pull the category label from the last ':' delimited part of column name
    df['category'] = [x.split(':')[-1] for x in full_name_list]
    # assign color digits according to category label
    cat_list = sorted(df['category'].unique())
    num_cats = len(cat_list)
    dic = dict()
    # use the presence of decimial point '.' to see if there's doses or not
    # assumed to be in second to last position if so
    if any(['.' in x for x in full_name_list]):
        hdrs = ['category', 'full', 'x', 'y', 'dose', 'color']
        df['dose'] = [x.split(':')[-2] for x in full_name_list]
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
    fig, ax = plt.subplots()
    for name, group in df.groupby('category'):
        ax.scatter(group['x'], group['y'], c=group['color'], label=name)
    ax.set_title(title.replace('_', ' ').replace('=',':'))
    #if len(df) < 16:
    #    ax.legend()
    plt.savefig(outpath + 'figs/' + title + '.png')
    for l, x, y in zip(df.index.values, df['x'].values, df['y'].values):
        plt.annotate(l, xy=(x, y), xytext=(-2,2),
                     textcoords='offset points', ha='right', va='bottom')
    plt.savefig(outpath + 'figs_annot/' + title + '.png')
    plt.close()
    df_save = df[['category', 'full', 'x', 'y']]
    df_save.to_excel(outpath + 'coords/' + title + '.xlsx')


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
        rdist[n] = math.sqrt((p1[0] - p2[0])*(p1[0] - p2[0]) + (p1[1] - p1[1])*(p1[1] - p1[1]))
    vals = list(dist.values())
    avg = round_num(np.mean(vals))
    med = round_num(np.median(vals))
    x_range = max(f['x'])
    y_range = max(f['y'])


    return avg, med



def main():
    path = '/Users/wrb/Dropbox/Areas of Focus/_Genometry/Analysis Projects/tSNE testing/'

    outpath = path + 'output/'

    try:
        os.makedirs(outpath)
        os.makedirs(outpath + 'coords')
        os.makedirs(outpath + 'figs')
        os.makedirs(outpath + 'figs_annot')
    except:
        pass

    # assumes the data file to be samples as columns, features as rows
    # single header row with plate:well:(batch):(dose):name/category eg NVA107:C02:A:0.0001:NVP001
    # the batch will just be present in full name where dose will convert to plotting diff sizes

    dat_path = '/Users/wrb/Dropbox/Areas of Focus/_Genometry/Analysis Projects/tSNE testing/gjb all 1+3um cons minus 2 loud.xlsx'

    print('loading ', dat_path)

    ld = {'ciclopirox': 'blue',
          'sirolimus': 'green',
          'tanespimycin': 'purple',
          'trichostatin': 'red',
          'vorinostat': 'orange',
          'test': 'grey',
          'null': 'black'}

    d = pd.read_excel(dat_path)

    # compose the target/label list
    inl = [x.split(':')[-1] for x in d.columns.values]
    # translate target/labels into colors per ld label dict defined above
    inl = ['test' if 'NVP' in x else x for x in inl]

    # temp switch to full column headers
    #labels = [ld[x] for x in inl]

    labels = d.columns
    #labels = [x.split(':')[-1][-4:] for x in d.columns.values]

    #print(labels)

    ds = scale(d)


    hdr = 'full'
    d = d.T
    run_methods(d, outpath, labels, hdr)
    setup_tsne(d, labels, outpath, hdr)

    hdr = 'scale'
    d = ds.T
    run_methods(d, outpath, labels, hdr)
    setup_tsne(d, labels, outpath, hdr)


if __name__ == '__main__':
    main()



#!/usr/bin/python -tt

"""
pa = plate analysis, a wide ranging collection of analytical functions for
investigating and manipulating datasets 

"""


import pandas as pd
import numpy as np
import gct, sys, os, csv, math, glob, pickle, gt, shutil, statistics
import collections as cll
import matplotlib.pyplot as plt
import seaborn as sns
import skyline, plottools, imgs, decimal
import dim_reduct
from sklearn import preprocessing
from sklearn.metrics.pairwise import euclidean_distances


def gene_shape_search(df, g, splot=True, top=50, pplot=False, dplot=False, ret=False):
    """ search for and plot genes closest to designated gene by scaled euclidean distance beween
    gene vectors. splot is summary plot of dist distribution, top limits results to that number of genes,
    pplot is to generate plate plats, gplots to generate dot plots """
    try:
        ptitle = df.name
    except:
        ptitle = df.index[0].split(':')[0]
    ds = scale_dataset(df)
    dist = comp_euclidean_dist(ds, ds.loc[g])
    if splot is True:
        fig, ax = plt.subplots()
        ax.set_title(ptitle + ' euclid distance from ' + g)
        ax.plot(dist.values[:500])
        ax.set_xlabel('genes')
        ax.set_ylabel('euclidean distance')
    if pplot is True or dplot is True:
        for i, g in enumerate(dist.index[:top]):
            if pplot is True:
                plottools.make_plateplot(df.loc[g], name=str(i) + '-' + ' plate')
                plt.close()
            if dplot is True:
                plottools.plot_gene(df.loc[g], name=str(i) + '-' + ' dot - ' + g)
                plt.close()
    if ret is True:
        return dist


def scale_dataset(df):
    """ uses scikit learn package to set each feature to have mean zero and standard deviation of 1 """
    res = preprocessing.scale(df)
    try:
        result = pd.DataFrame(data=res, index=df.index, columns=df.columns)
    except:
        print('trying for series')
        result = pd.DataFrame(data=res, index=df.index)
    return result


def comp_euclidean_dist(df, vect, top=False):
    """ pass in a gene vector and dataframe to compare pairwise euclidean distances with the vector
    and return a ranked list of the top n number of closest features """
    try:
        dat = euclidean_distances(df, vect.values.T)
    except:
        dat = euclidean_distances(df, vect.values[np.newaxis])
    #print(dat.flatten())
    dist = pd.Series(data=dat.flatten(), index=df.index)
    dist.sort_values(inplace=True)
    if top is False:
        return dist
    else:
        return dist.iloc[:top]


def predict_cells(input, save=False):
    """ can accept directory and loop through files or one dataframe at a time,
     uses v1.0 of the SVM classifier to consolidate reps to consensus and return prediction
     when save is True a dataframe will be saveh"""

    with open('/Users/WRB/Dropbox/bin/python/celllineclassifier.p', 'rb') as file:
        clf = pickle.load(file)
    if isinstance(input, str):
        if os.path.isdir(input):
            vlist = gt.globit(input, '*_Qctrl_n*')
            if len(vlist) == 0:
                vlist = gt.globit(input, '*QNORM*')
        else:
            vlist = [input]
    elif isinstance(input, pd.Series):
        try:
            res = clf.predict([input])[0]
        except:
            print('error with series prediction')
            res = None
        return res
    else:
        vlist = input
    res_table = pd.DataFrame()
    for f in vlist:
        try:
            d, h = gct.extractgct(f)
        except:
            vlist[0] = d
            vlist[1] = h
        ds, hs = gt.dsub(d, h, {'type':'vehicle'})
        if len(ds) == 0:
            print('error, maybe using ZS file? use QNORM instead')
            return None
        for b in hs.batch.unique():
            dsb, hsb = gt.dsub(ds, hs, {'batch':b})
            med = dsb.median(axis=1).values
            shn = gt.get_shn(f) + '-' + b
            res = clf.predict([med])[0]
            res_table.loc[shn,'cell'] = res
            print(f'{shn} - {res}')
    if save is True:
        res_table.to_csv(gt.dflt_outpath(fn='cell_predictions.csv'), sep='\t')
    return res_table


def plot_veh_matrices(path='dflt', batch='A', up=.75, low=0, kind='euclid', getcells=True):
    """ gather data for vehicle matrix and plot it regular and clustered in euclidean """

    d, h = get_vehicle_matrix(path=path, batch=batch, getcells=getcells)

    h = gt.gen_label(h, 'pb')

    dim_reduct.make_dendrogram(d, labels=h.label)

    plottools.plot_correlation_matrix(d, h, lower=low, upper=up, ptype=kind,
                                      title=f'vehicles {kind}', outpath=True, sparselabel=False)

    plottools.plot_clustered_correl(d, h, ptype=kind,
                                    title=f'clustered vehicles {kind}', outpath=True)


# run this on the overall folder, looking for Qtcrl or QNORM
def get_vehicle_matrix(path='dflt', batch='all', delim=':', getcells=False):
    """" for the path load all files and collapse vehicles, plot matrix
    batches can be all or 'A' only to just take the first one. getcells will re-predict cells """
    path = gt.check_dfltarg(path, os.path.join(gt.check_desktop(), 'newQC'))

    flv = gt.globit(path, '*Qctrl*')
    if len(flv) == 0:
        flv = gt.globit(path, '*_QNORM_*')

    # should put in a check to extract from regular qnorms
    dlist, hlist = [], []
    for f in flv:
        d, h = gct.extractgct(f)
        h['plate'] = h['plate'].apply(lambda x: x[:6])
        d, h = gt.dsub(d, h, {'type':'vehicle'})
        if batch == 'all':
            for b in h.batch.unique():
                ds, hs = gt.dsub(d, h, {'batch':b})
                med = ds.median(axis=1)
                hs = gt.gen_label(hs, 'pb', delim=delim)
                dlist.append(med)
                hlist.append(hs.iloc[0])
        elif batch == 'A':
            ds, hs = gt.dsub(d, h, {'batch': 'A'})
            med = ds.median(axis=1)
            hs = gt.gen_label(hs, 'pb', delim=delim)
            dlist.append(med)
            hlist.append(hs.iloc[0])
        else:
            med = d.median(axis=1)
            hs = gt.gen_label(hs, 'p', delim=delim)
            dlist.append(med)
            hlist.append(hs.iloc[0])

    vdf = pd.concat(dlist, axis=1)
    vh = pd.DataFrame(hlist)
    vdf.columns = vh.label
    if getcells is True:
        vh['cell2'] = vh.label.apply(lambda x: predict_cells(vdf[x]))
        vh['label'] = vh.label + delim + vh.cell2
    vdf.columns = vh.label
    return vdf, vh


# currently plotting euclidean, not straight up correlation between plates
def plate_comparison(flist, scat=True, corr=True, dat=False):
    """ for a range of gct files consolidate median, mean, std, cv from each and then
     create pairwise scatterplots for each. flexible in number of plates """
    if isinstance(flist,str):
        flist = gt.splitpaths(flist, '.gct')
    outpath = gt.dflt_outpath(fldr_name='comparisons')
    ddict = cll.OrderedDict()
    hdict = cll.OrderedDict()
    for i,f in enumerate(flist):
        name = gt.get_shn(f)
        df, h = gct.extractgct(f)
        ddict[name], hdict[name] = df, h
        if i == 0:
            baseindex = df.index
    medians = pd.DataFrame(index=baseindex)
    medians.name = 'median gene values'
    stdev = pd.DataFrame(index=baseindex)
    stdev.name = 'gene standard deviations'
    cv = pd.DataFrame(index=baseindex)
    cv.name = 'gene coefficient of variation'
    average = pd.DataFrame(index=baseindex)
    average.name = 'gene average'
    for n, d in ddict.items():
        medians[n] = d.median(axis=1)
        stdev[n] = d.std(axis=1)
        cv[n] = d.std(axis=1) / d.mean(axis=1)
        average[n] = d.mean(axis=1)
    for dset in [medians, stdev, cv, average]:
        if scat is True:
            sns.pairplot(dset)
            plt.tight_layout()
            plt.suptitle(dset.name)
            plt.savefig(os.path.join(outpath, dset.name+'scatter.png'))
            plt.close()
        if dat is True:
            dset.to_excel(os.path.join(outpath, dset.name + '.xlsx'))
        if corr is True:
            ax = plottools.plot_euclidean(dset, dset.columns)
            ax.set_title(dset.name)
            plt.tight_layout()
            plt.savefig(os.path.join(outpath, dset.name + 'matrix.png'))
            plt.close()


def plot_lmconcs(flist, stack=False, test=False):
    """ plot lm concs with full output for all files in list, w/ optional joins"""
    if isinstance(flist, str):
        flist = gt.splitpaths(flist,'.gct')
    for f in flist:
        d, h = gct.extractgct(f)
        #ds, hs = gt.dsub(d, h, {'name':['5-Iodotubercidin', 'ERK5-IN-1']})
        outpath = gt.dflt_outpath(fldr_name='landmark concs')
        plottools.plot_landmark_concs(d, h, genes='all', labels='wells', outpath=outpath, test=test)
    if stack is True:
        imgs.bulk_stack(outpath, orient='vert', delim='_', idx=2)


def plot_skylines(df, h, argdict, title='dflt'):
    """ pass a dataframe and argument dictionary to just plot skylines matching criteria
     if a dict key is passed with True as the arguments, it will plot all conditions from that category
     eg {'name':myname, 'dose':True } will plot all doses of that name """

    cats = [x[0] for x in argdict.keys()]
    cats = ''.join(cats)
    print(cats)

    myargdict = {}
    for k, v in argdict.items():
        if isinstance(v, bool) and v is True:
            continue
        myargdict[k] = v

    ds, hs = gt.dsub(df, h, myargdict)

    ddict = breakdown(ds, hs, cats, dic=True)

    if isinstance(title, str) and title == 'dflt':
        title = df.name

    for name, dat in ddict.items():
        if isinstance(dat, bool) and dat is True:
            continue
        mytitle = title + ' - ' + name
        skyline.new_skyline(dat, title=mytitle)


def assemble_consensus(df, h, cats, sc=True, legend='brief', plot=False, skyl=False, n=None, save=False,
                       ret=True, test=False):
    """ tool to assemble replicate zscore consensus, pass df, header and the breakdown categories 'nd' for instance
    will return the consolidated df and header file

    can pass in multiple gct files as a single string via F6 shortcut

    ccs will calculate the zscore correlation of replicates, and insert that into header df
    plot will use seaborn pairplot to visualize the calculated rep correlations above
    skyl controls skyline plot generation, can be True to plot all ind reps plus consensus
    n argument is a limiter to only consider treatments with enough replicates, including into consensus gct!!
    save will save the consensus gct file
    """
    
    if isinstance(df, str):
        df, h = gct.extractgct(df)
    else:
        print('error in loading dataframe')

    outpath = gt.dflt_outpath(fldr_name='output figs')
    try:
        pname = df.name
    except:
        pname = h.addr[0].split(':')[0]
    outpath = os.path.join(outpath, pname)
    try:
        os.mkdir(outpath)
    except:
        pass

    subs = breakdown(df, h, cats, dic=False)

    con_data = pd.DataFrame(index=df.index)
    addnl = []
    addnl.extend(['corr', 'all ccs'])
    addnl.extend(['prom', 'all proms', 'porder'])

    if addnl != []:
        con_header = pd.DataFrame(index=list(h.columns.values)+addnl)
    else:
        con_header = pd.DataFrame(index=h.columns)

    for ds, hs in subs:
        if n is not None:
            if len(ds.columns) < n:
                print('not enough reps', hs.iloc[0])
                continue

        c = consensus(ds, name='first')
        con_data = pd.concat([con_data, c], axis=1)

        new_annot = hs.iloc[0,:].copy().T
        new_annot.well = hs['well'].values
        new_annot.addr = hs['addr'].values


        corrs = []
        for i in range(len(ds.columns)):
            for j in range(1 + i, len(ds.columns)):
                corrs.append(round(ds.iloc[:, i].corr(ds.iloc[:, j], method='pearson'), 2))
        if len(corrs) == 0:
            # print('corrs = na')
            # print(hs.iloc[0].values)
            new_annot['corr'] = np.nan
            new_annot['all ccs'] = np.nan
        elif len(corrs) == 1:
            new_annot['corr'] = round(corrs[0],2)
            new_annot['all ccs'] = corrs
        else:
            new_annot['corr'] = round(np.percentile(corrs, 75), 2)
            new_annot['all ccs'] = corrs
        corrs = [decimal.Decimal(x) for x in corrs]
        new_annot['corr'] = pd.to_numeric(new_annot['corr'])

        proms = abs(ds).sum(axis=0).round().values
        porder = hs['well'].values
        new_annot['prom'] = round(np.percentile(proms, 75))
        new_annot['all proms'] = proms
        new_annot['porder'] = porder

        if plot is True:
            ds.columns = [x + ' - ' + hs.loc[x]['batch'] for x in ds.columns]
            ax = sns.pairplot(ds)
            myoutpath = os.path.join(outpath, 'rep zs scatter')
            try:
                os.mkdir(myoutpath)
            except:
                pass
            plt.savefig(os.path.join(myoutpath, h.plate[0] + '-' + ds.name + '.png'))
            plt.close()

        con_header = pd.concat([con_header, new_annot], axis=1)

        if skyl is True:
            myoutpath = os.path.join(outpath, 'skyline')
            try:
                os.mkdir(myoutpath)
            except:
                pass
            try:
                name = hs.iloc[0]['name'] + '-' + str(hs.iloc[0]['dose']) + '-' + hs.iloc[0]['batch']
            except:
                name = hs.iloc[0]['name'] + '-' + hs.iloc[0]['batch']
            name = name.replace('.', ',')
            title = pname + '-' + name
            myoutpath = os.path.join(myoutpath, title)
            skyline.new_skyline(ds, title=title, outpath=myoutpath)

        if test is True:
            break

    con_header = con_header.T

    if sc is True:
        try:
            pname = df.name
        except:
            pname = h.addr[0].split(':')[0]
        title = pname + ' sc plot'
        outpath = gt.dflt_outpath(fn=pname+'_scplot.png')
        kwargs = {'x': 'corr', 'y': 'prom',  'data': con_header}

        kwargs.update({'alpha': .75, 'style':'type', 'legend': legend})

        if 'd' in cats:
            kwargs['hue'] = 'name'
            kwargs['size'] = 'dose'
            kwargs['sizes'] = (40, 400)

        # this is experimental
        else:
            kwargs['sizes'] = (50)
            kwargs['hue'] = 'name'

        g = sns.relplot(**kwargs)
        g.fig.suptitle(title)
        g.fig.set_size_inches(7,5)
        if legend is not None:
            for lh in g._legend.legendHandles:
                lh.set_alpha(.75)
        g.savefig(outpath, bbox_inches='tight')
        plt.close()

        con_header = gt.gen_label(con_header, 'nb')
        newfig, newax = dim_reduct.seaborn_scatter(con_header, title, outpath, x='corr', y='prom', ptype='ax', save=False)
        dim_reduct.html_scatter(newfig, newax, con_header, 'corr', 'prom', title)
        plt.close()

    con_data.name = df.name

    if save is True:
        gct.save_headergct(con_data, con_header, gt.dflt_outpath(fn=df.name+'_consensus.gct'))
    if ret is True:
        return con_data, con_header


def consensus(ds, name='dflt'):
    """ merges ZScore instances in a passed dataframe into a conservative consensus - min abs value
    as long as all reps agree on direction, otherwise zero. uses a fancy np.select method to set values
    if name is 'dflt' the column name will be: FPA001:K03-N04_(3)
    otherwise if name is 'first' name will just be first well FPA001:K03 (better for header integration) """
    choices = [ds.min(axis=1), ds.max(axis=1)]
    conds = [(ds > 0).all(axis=1), (ds < 0).all(axis=1)]
    newdf = pd.Series(data=np.select(conds, choices, default=0), index=ds.index)
    if name is 'dflt':
        newdf.name = ds.columns[0] + '-' + ds.columns[-1].split(':')[-1] + '_(' + str(len(ds.columns)) + ')'
    elif name is 'first':
        newdf.name = ds.columns[0]
    return newdf


def breakdown(df, h, cats, dic=True, genes=None):
    """ takes a dataframe and header and the categories to break down by 'b' batch, 'c' cell, 'n' name, 'd' dose.
    returns a dictionary with the key as the description and the dataframe as the value.
    'w' is also supported as breakdown by well - useful for many plates with identical layout

    if dic is True a dictionary is returned, with a key title and dataframe value
    if dic is False then list is returned, of tuples with dataframe and header

    """

    if genes is not None:
        genes = gt.get_genes(genes)
        df = df.loc[genes]

    if 'd' in cats:
        try:
            dose_col = [x for x in h.columns if 'dose' in x or 'dilution' in x][0]
        except IndexError:
            print('dose column error')
    else:
        dose_col = None

    vd = cll.OrderedDict()
    subs = []

    cd = {'c': 'cell',
          'b': 'batch',
          'd': dose_col,
          'n': 'name',
          'w': 'well',
          'p': 'plate'}

    clist = []

    for c in cats:
        try:
            clist.append(cd[c])
        except IndexError:
            print('error, more than 3 categories')

    cat1 = clist[0]
    group1 = sorted(h[cat1].dropna().unique())
    for e1 in group1:
        argdict = {cat1: e1}
        try:
            cat2 = clist[1]
            for e2 in sorted(gt.hsub(h, {cat1: e1})[cat2].dropna().unique()):
                argdict.update({cat2: e2})
                try:
                    cat3 = clist[2]
                    for e3 in sorted(gt.hsub(h, {cat1: e1, cat2: e2})[cat3].dropna().unique()):
                        argdict.update({cat3: e3})
                        hdr = f'{e1}-{e2}-{e3}'
                        if dic is True:
                            vd.update({hdr: gt.dosub(df,h, argdict, name=hdr)})
                        else:
                            subs.append(gt.dsub(df,h, argdict, name=hdr))
                except IndexError:
                    hdr = f'{e1}-{e2}'
                    if dic is True:
                        vd.update({hdr: gt.dosub(df, h, argdict, name=hdr)})
                    else:
                        subs.append(gt.dsub(df, h, argdict, name=hdr))
        except IndexError:
            hdr = f'{e1}'
            if dic is True:
                vd.update({hdr: gt.dosub(df, h, argdict, name=hdr)})
            else:
                subs.append(gt.dsub(df, h, argdict, name=hdr))

    if dic is True:
        return vd
    else:
        return subs


def rank_files(path):
    """ reads a zscore gct file (probably consensus) and saves file of ranks
            of each gene within each sample - from highest to lowest values"""
    flist = gt.get_flist(path, '.gct')
    for f in flist:
        g = gct.Gct(f)
        print(f)
        d = g.build_dframe()
        ds = d.rank(ascending=False)
        outpath = os.path.join(os.path.split(path)[0], g.shortname + '-ranks.gct')
        gct.save_headergct(ds, outpath, g.fh)


def grab_ranks(path, feat, hilo=1, t=40):
    """ survey folder for the wells in whic given gene istop ranked wells of given feat by sorted z score). generates
        overall list ofdefault rank output is in descending order, highest zscore = 1
        hilo: 1 = high, upregulated genes (default rank order)
                    0 = low, downnregulated genes """
    outpath = os.path.join(path, '_rank_summary.txt')
    flist = gt.get_flist(path, 'ranks.gct')
    # set dummy starting point for low rank
    lowest = 500
    # create blank template dataframe
    summary = pd.DataFrame()
    for f in flist:
        d, h = gct.extractgct(f)
        # flip rank order as needed
        if hilo > 1:
            d = 978 - d
        # get column ids for ranks below threshold
        wells = d.columns[d.ix[feat] < t]
        # extract portion of dataframe
        ranks = d.ix[feat, wells]
        ranks = pd.DataFrame(ranks)
        # assign plate column to each well id entry, re-order cols
        ranks['plate'] = gt.get_shn(f).split('-')[0]
        # concat portion to overall dataframe
        summary = pd.concat([summary, ranks])
        # check and store the lowest rank
        newlow = min(d.ix[feat])
        if newlow < lowest:
            lowest = newlow
    # re-shuffle the column order
    summary['well'] = summary.index.values
    summary = summary[['plate', 'well', feat]]
    print('\n', feat, int(lowest))
    summary.to_csv(outpath, sep='\t', index=None)


def get_zscore(fpath, save=True, my_mad=None):
    """ merged from separate zscore file. can either save the resulting file or return data
    the first fpath argument can be a file path or a [d, h] object already"""
    # basic setup
    if isinstance(fpath, str):
        g = gct.Gct(fpath)
        g.get_headers()
        df, h = gct.extractgct(fpath)
    else:
        try:
            df = fpath[0]
            h = fpath[1]
        except:
            print('error with path')

    zsd = cll.defaultdict(dict)
    pname = gt.get_shn(fpath)

    for b in h['batch'].dropna().unique():
        if b == 'na':
          continue
        print('running zscore for {} batch {}'.format(pname, b))
        vw = gt.hsub(h, {'batch':b, 'type':'vehicle'}).index.values
        if len(vw) == 0:
            break
        veh = df[vw]
        # get median value across vehicle populations
        med = veh.median(axis=1)

        # populate the absolute deviation values per gene
        ad = cll.defaultdict(list)
        for v in veh.columns:
            for f in veh.index:
                ad[f].append(abs(med[f] - veh[v][f]))
        # assemble the median absolute value per gene
        mad = {}
        for k, v in ad.items():
            r = statistics.median(v)
            if 0 < r < 0.1:
                r = 0.1
            mad[k] = r
        # using the above progress though test and poscon wells
        # to calculate sample zscores
        tw = list(h[(h['batch'] == b) & (h['type'] == 'test')].index.values)
        pw = list(h[(h['batch'] == b) & (h['type'] == 'poscon')].index.values)
        wells = tw + pw
        for w in df[wells].columns:
            for feat in df.index:
                if my_mad is not None and mad[feat] < my_mad:
                    zs = (df[w][feat] - med[feat]) / (my_mad * 1.486)
                elif mad[feat] == 0:
                    zs = 0
                else:
                    zs = (df[w][feat] - med[feat]) / (mad[feat] * 1.486)
                zsd[w][feat] = '{0:.3f}'.format(zs)

    # transform into dataframe, set index, null nonsense
    zsdf = pd.DataFrame(zsd)
    hs = h.loc[zsdf.columns]
    zsdf = zsdf.replace(['inf', '-inf'], np.nan).fillna('nan')
    if save is True:
        outpath = '{}_ZS.gct'.format(fpath.split('_', 1)[0])
        gct.save_headergct(zsdf, hs, outpath)
    else:
        return zsdf, hs


def get_dir_zscores(opath):
    # test whether input is file or directory, do single or loop through
    # all files as appropriate
    if os.path.isfile(opath) and opath.endswith('.gct'):
        get_zscore(opath)
    elif os.path.isdir(opath):
        for file in os.listdir(opath):
            if file[0] != '.':
                if file.endswith('.gct'):
                    targfile = os.path.join(opath, file)
                    get_zscore(targfile)



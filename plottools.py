#!/usr/bin/python

""" repository for subfunctions and plotting functions to generate plots of genes,
plates etc, all available as standalone functions """


import gct, os, gt, string, scipy, math, sys, pa, dim_reduct
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import collections as cll
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from matplotlib import rcParams


def plot_inv_level(flist, inv_level, scale=False):
    """ create plate plot for a given invariant level across plates, to see if any pattern
    can be either fullqnorm or reg qnorm """
    flist = gt.splitpaths(flist, ext='.csv')
    for f in flist:
        d = gcsv.open_as_gct(f)
        shn = '_'.join(os.path.split(f)[-1].split('_')[:2])
        if scale is False:
            vctr = d.loc[f'0_INV_{inv_level}']
        elif scale is True:
            vctr = d.loc[f'0_INV_{inv_level}'] / d.loc['0_INV_10']
        make_plateplot(vctr, name=shn + '-lvl' + str(inv_level))



def arrange_cohorts(d, h, genes, dtype='auto', title='auto', n=False, size=25):
    """ takes dataframe and genes and plo"""
    if title == 'auto':
        title = h.index[0].split[':'][0]
    if n is not False:
        for g in genes:
            for n in h['name'].unique():
                vdict = {}
                for c in h['cell'].unique():
                    try:
                        vdict[c] = gt.dosub(d[h.index], h, {'name': n, 'cell': c}).loc[g]
                    except:
                        vdict[c] = [0]
                mytitle = title + ' ' + n + ' ' + g
                print(mytitle)
                plottools.plot_cohorts(vdict, dtype=dtype, title=mytitle, incr=5, size=size)
    elif n is False:
        for g in genes:
            vdict = {}
            for c in h['cell'].unique():
                vdict[c] = gt.dosub(d[h.index], h, {'cell': c}).loc[g]
            mytitle = title + ' ' + g
            print(mytitle)
            plottools.plot_cohorts(vdict, dtype=dtype, title=mytitle, incr=5, size=size)



def compare_plate_genes(flist, genelist, numrows=3, type=True, plate=False, title='dflt', outpath='dflt', remove=True):
    """ plots the listed genes across the dataframes provided in list of file paths flist.
     the plots will be generated and then combined using img bulk stack. orient is direction
     of the joined images. type will include sample type color coding -- should add in grid support ---  """
    if isinstance(flist,str):
        flist = gt.splitpaths(flist, '.gct')
    if outpath == 'dflt':
        outpath = gt.dflt_outpath(fldr_name='tmp_imgs')
    for f in flist:
        d, h = gct.extractgct(f)
        for g in genelist:
            if plate is False:
                plot_gene_wtypes(d.loc[g], h, name=d.name+'-'+g, outpath=outpath)
            elif plate is True:
                make_plateplot(d.loc[g], name=d.name + '-' + g, outpath=outpath)
    if title != 'dflt':
        combined_outpath = gt.dflt_outpath(fldr_name=title + ' combined_imgs')
    else:
        combined_outpath = gt.dflt_outpath(fldr_name='combined_imgs')
    if numrows is False:
        imgs.bulk_stack(outpath, outpath=combined_outpath, delim='-', idx=1, pad=.05)
    else:
        imgs.bulk_stack_grid(outpath, outpath=combined_outpath, numrows=numrows, delim='-', idx=1, pad=.05)
    if remove is True:
        shutil.rmtree(outpath)

def blue_red_cmap(colors=['cornflowerblue', 'white', 'red']):
    """ define continuous color gradation from passed colors to be used for plotting
    royalblue is a nice deeper blue, cornflower blue is lighter """

    cmap = LinearSegmentedColormap.from_list('mycmap', colors)

    return cmap


def sweep_matrices(df, h, title='dflt', kind='reg'):
    """ higher level function to create range of trimmed matrices in pearson correlation and euclidean
    space, default type is just regular unclustered matrix """

    uppers = [1, .9, .75, .6]
    lowers = [0, 0.1, .25, .35]

    ptype = ['euclid', 'corr']

    fname = df.columns[0].split(':')[0]

    for pt in ptype:
        for u in uppers:
            for l in lowers:
                title = f'{fname}-{pt}-{u}{l}'
                plot_correlation_matrix(df, h, title=title, ptype=pt, lower=l, upper=u, outpath=True)



def plot_clustered_correl(df, h, title='dflt', ptype='euclid', outpath='dflt', lower=0, upper=0.75):
    """ use ward's minimum distance linkage from sci py dendrogram to determine sample ordering
    for plotting in a pairwise distance matrix """

    try:
        h = h.set_index('label')
        h['label'] = h.index
        df.columns = h.label
    except:
        pass

    corder = dim_reduct.clustered_order(df)
    ord_dict = dict(zip(range(0,len(df.columns)), df.columns))

    sample_order = list([ord_dict[x] for x in corder])

    newdf = df[sample_order]

    h = h.reindex(sample_order)

    plot_correlation_matrix(newdf, h, ptype=ptype, title=title, upper=upper, lower=lower,
                            sort=False, sparselabel=False, outpath=True)


def plot_correlation_matrix(df, h, ptype='corr', title='dflt', labels='dflt', sort=False, lower=0.25,
                            upper=1.0, outpath=False, cmap='dflt', sparselabel=False, grid=False):
    """ plots pearson correlation matrix between columns of the passed in dataframe. the labels can be used
     to sort the samples, sparselabel only prints one category label per section/cluster, and outpath will save
     to designated location, otherwise just display
     the 'lower' argument trims bottom of graph, so that there's less noise at the bottom end

     the sparselabel designation requires sorting, otherwise things work out funny
     _should improve label handling """

    fig, ax = plt.subplots()
    fig.set_size_inches(8.5, 8.5)

    if isinstance(labels, str) and labels == 'dflt':
        try:
            mylabels = h.label.values
        except:
            mylabels = df.columns
    if labels is None:
        mylabels = len(df.columns) * ['']
    else:
        mylabels = labels

    if sparselabel is True:
        sort = True

    if sort is True:
        print(f'{len(df.columns)} columns, {len(labels)} labels')
        keyd = dict(zip(df.columns, mylabels))
        neword = sorted(df.columns, key=lambda x: keyd[x])
        mylabels = [keyd[x] for x in neword]
        df = df[neword]

    if ptype == 'corr':
        corr = df.corr()
    elif ptype == 'euclid':
        cmap = 'rev'
        corr = get_euclidean(df, df=True)
        max = corr.max().max()
        upper = (1-lower) * max
        lower = 0

    if lower is not None:
        corr = corr.clip(lower=lower)
    if upper is not None:
        corr = corr.clip(upper=upper)

    if cmap == 'dflt':
        cmap = blue_red_cmap()
    elif 'rev' in cmap:
        cmap = blue_red_cmap(['red','white','cornflowerblue'])


    cax = ax.imshow(corr, interpolation='nearest', cmap=cmap)
    #cbar = fig.colorbar(cax, ticks=[-1,0,1])

    if sparselabel is False:
        # x axis
        minor = np.arange(0.5, len(df.columns), 1)
        major = np.arange(0.5, len(df.columns), 1)
        ax.set_xticks(major, minor=False)
        ax.xaxis.set_tick_params(size=0)
        ax.set_xticks(minor, minor=True)
        # y axis
        minor = np.arange(0.5, len(df.columns), 1)
        major = np.arange(0, len(df.columns), 1)
        ax.set_yticks(major, minor=False)
        ax.yaxis.set_tick_params(size=0)
        ax.set_yticks(minor, minor=True)

        ax.set_xticklabels(mylabels, rotation=45, ha='right')
        ax.set_yticklabels(mylabels)

    elif sparselabel is True:
        unq_labels = sorted(list(set(mylabels)))
        cntr = cll.Counter(mylabels)
        ticks, mylabels, i, major, minor = [1], [], -.5, [], []
        for cat in unq_labels:
            chunk_size = cntr[cat]
            label_loc1 = i + chunk_size/2
            label_loc2 = i + chunk_size
            #print(f'{cat} : len {cntr[cat]} at position {label_loc1}')
            major.append(label_loc1)
            minor.append(label_loc2)
            #ticks.extend([label_loc1, label_loc2])
            #mylabels.extend(['', cat])
            mylabels.append(cat)
            i += chunk_size

        ax.set_xticks(major, minor=False)
        ax.xaxis.set_tick_params(size=0)
        ax.set_xticks(minor, minor=True)
        ax.set_yticks(major, minor=False)
        ax.yaxis.set_tick_params(size=0)
        ax.set_yticks(minor, minor=True)

        ax.set_xticklabels(mylabels, rotation=45, ha='right')
        ax.set_yticklabels(mylabels)

        if grid is True:
            ax.grid(which='minor', axis='both', color='black')

    if title == 'dflt':
        try:
            title = df.name
        except AttributeError:
            title = df.columns[0].split(':')[0]

    ax.set_title(f'{title} - n={len(df.columns)} corr matrix', style='oblique')

    plt.tight_layout()

    if outpath is not False:
        if outpath is True:
            plt.savefig(gt.dflt_outpath(fn=title + '_corr.png'))
        else:
            if '.png' in outpath:
                plt.savefig(outpath)
            else:
                plt.savefig(outpath + '.png')
        plt.close()


def plot_landmark_concs(df, h, maxy=12, cats='n', labels='dflt', genes='test100', outpath='dflt',
                        title='dflt', dosenum='dflt', test=False):
    """ plot many or all landmarks, should pass in a subset dataframe and header which
    should be the consensus ZS file. can contain many different names + doses, will auto breakdown by 'nd'
    a single line per gene is plotted for the ZS across all concentrations
     labels can be 'dflt' for just incr numbers, or 'wells' for address, or 'dose' for numbers """
    # txt_args = {'fontsize': 8,
    #             'rotation': 90,
    #             'fontweight': 'bold'}

    if outpath is 'dflt':
        outpath = gt.dflt_outpath()
    df, h = gt.dsub(df, h, {'type':'test'})
    names = h.name.dropna().unique()
    doses = gt.hsub(h, {'name': names[0]})['dose'].dropna().unique()
    if len(gt.hsub(h, {'name': names[0], 'dose': doses[0]})) > 1:
        print('dataframe not collapsed to consensus, bogus lm concs')
        print(gt.hsub(h, {'name': names[1], 'dose': doses[0]}).head())
    for ds, hs in pa.breakdown(df, h, cats, dic=False):
        #hs['dose'] = pd.to_numeric(hs['dose'])
        hs.sort_values('dose', ascending=True, inplace=True)
        ds = ds[hs.index]
        xrange = len(hs.dose.unique())-2
        ax = format_concentration_plot(xrange, maxy=maxy, width=4)
        ax.tick_params(axis='x', bottom='on', top='off', labelbottom='on')
        if dosenum == 'dflt':
            dose_range = range(len(hs.dose.unique()))
        else:
            dose_range = range(dosenum)
        ax.set_xticks(dose_range)
        if labels == 'dflt':
            ax.set_xticklabels([str(x + 1) for x in dose_range])
        elif labels == 'wells':
            # temporary labels
            ax.set_xticklabels(hs.index, rotation=45)
        elif labels == 'dose':
            ax.set_xticklabels(hs['dose'].unique(), rotation=45)
        else:
            try:
                ax.set_xticklabels(labels)
            except:
                print('problem with x range labels')

        # set title and name
        if title == 'dflt':
            try:
                mytitle = df.name
            except:
                mytitle = hs['plate'].values[0]
        mytitle = mytitle.strip('_sub')
        suffix = ''
        for c in cats:
            cat = gt.cats_lookup(c)
            attr = hs[cat].values[0][0]
            suffix += f' - {attr}'
        mytitle += suffix

        ax.set_title(mytitle, fontsize=14)
        for g in gt.get_genes(genes, df=df):
            data = ds.loc[g,:]
            ax.plot(data.values, linewidth=0.3)
        plt.tight_layout()
        plt.savefig(os.path.join(outpath, mytitle + '.png'))
        plt.close()
        if test is True:
            print('stopping after one iteration')
            break


def get_euclidean(data, df=False):
    """ transforms a dataframe into a pairwise euclidean distance matrix """
    Y = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(data.T, 'euclidean'))
    if df is False:
        return Y
    else:
        return pd.DataFrame(Y, columns=data.columns, index=data.columns)


def plot_euclidean(df, labels='dflt', upper=None, rot=None, fontsize=None, tick_denom=1, dat=None):
    ''' takes dataframe and plots square euclidean distance, labels w/ labels
    must define ax = plot_euclidean() to get it back an operable object'''
    if dat is not None:
        Y = dat
    else:
        #Y = get_euclidean(df.T)
        Y = get_euclidean(df)

    if upper is not None:
        Y = Y.clip(max=upper)

    fig, ax = plt.subplots(figsize=(12,12))

    ax = sns.heatmap(Y, square=True, cmap='bwr_r', ax=ax)

    if tick_denom == 1:
        ticks = np.arange(0.5, len(df.columns), 1)
    else:
        ticks = np.arange(0, len(df.columns), tick_denom)

    if isinstance(labels, str):
        if labels == 'dflt':
            labels = df.columns
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    v_align = 'center'

    if rot is None:
        rot = 90
        hor_align = 'center'
    else:
        hor_align = 'center'

    if tick_denom != 1:
        hor_align = 'center'
        v_align = 'center'

    if fontsize == None:
        ax.set_xticklabels(labels, rotation=rot, ha=hor_align)
        ax.set_yticklabels(labels, rotation=0, va=v_align)
    else:
        ax.set_xticklabels(labels, rotation=rot, ha=hor_align, fontsize=fontsize)
        ax.set_yticklabels(labels, rotation=0, va=v_align, fontsize=fontsize)


    try:
        name = df.name + ' euclidean distance'
    except:
        name = df.columns[0].split(':')[0]
    ax.set_title(name)
    #plt.colorbar(ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    #plt.suptitle(name)
    return ax


def gen_euclideans(df, labels='dflt', rot=None, tick_denom=1, test=False):
    """ loops plotting euclidean matrix to use different upper trim boundaries and font sizes for labels """

    outdir = gt.dflt_outpath(fldr_name='matrices')

    try:
        name = df.name
    except:
        name = df.columns[0].split(':')[0]

    Y = get_euclidean(df)

    maxv = round(Y.max())

    lims = [1, .75, .5, .3, .15]

    for fs in [8, 5]:
        for ul in lims:
            cap = int(round(maxv * ul))
            ax = plot_euclidean(df, labels=labels, upper=cap, fontsize=fs, dat=Y, tick_denom=tick_denom, rot=rot)
            outpath = os.path.join(outdir, name + f'_euclidean_ul{str(ul).replace(".",",")}-fs{fs}.png')
            #plt.savefig(outpath, bbox_inches='tight')
            plt.savefig(outpath)
            plt.close()
            if test is True:
                sys.exit('test mode, quitting after one')


def fix_ticklabels(ax, l):
    ax.set_xticklabels(l, rotation=40, ha='right')
    #ax.set_yticklabels(l[::-1], rotation=0, va='center')
    ax.set_yticklabels(l, rotation=0, va='center')
    return ax


def set_titles(ax, title, xlabel, ylabel, xrot=0, yrot=0):
    """ shortcut to label axis with fig title, axis labels and set x+y tickmark rotation """
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=xrot)
    ax.set_yticklabels(ax.yaxis.get_majorticklabels(), rotation=yrot)
    

def format_concentration_plot(maxx, ptype='zs', width='dflt', height='dflt', maxy=10):
    """ basic simple y-axis only plot formatting. set the ylims externally
    (so this is flexible for zs and gct) """
    if width is 'dflt':
        width = maxx / 20
        if width > 14:
            width = 14
    if height is 'dflt':
        height = 3
    fig = plt.figure(figsize=[width, height])
    ax = fig.add_subplot(111)
    ax.set_facecolor("white")
    # format graph
    if ptype is 'zs':
        ax.set_ylim([-maxy, maxy])
        ax.set_yticks([-maxy, 0, maxy])
        #lo = '{:+>2}'.format(-maxy)
        #hi = '{:>2}'.format('{:+.0f}'.format(maxy))
        #ax.set_yticklabels([lo, 0, hi])
        #ax.set_yticklabels([-maxy, 0, '+' + str(maxy)])
    elif ptype is 'gct':
        ax.set_ylim([2, 16])
        ax.set_yticks(np.arange(2, 18, 2))
    elif ptype is 'csv':
        ax.set_ylim([2, 16])
        if maxy > 500:
            ax.set_yticks(np.arange(0, 10000, 2))
        else:
            ax.set_yticks(np.arange(0, 10000, 2))
    ax.set_xlim([0, maxx+1])

    # set the tick parameters
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(True)
    ax.tick_params(axis='x', bottom='off', top='off', labelbottom='off')
    ax.tick_params(axis='y', right='off', left='on')
    ax.tick_params(axis='y', direction='in', length=3, width=1)
    ax.yaxis.set_label_position('left')
    return ax


def find_plot_genes(df, thresh='dflt', lim=75):
    """ open a gct file and make list of genes which have ZS above the threshold, now updated to
     auto determine the list at 97.5% and above, also can input a hard value instead"""
    print(df.name)
    # auto determine threshold, currently 97.5
    thresh = gt.check_dfltarg(thresh, df.quantile(0.975).max())
    # apply filter
    subset = df[abs(df) > thresh]
    subset.dropna(axis=0, how='all', inplace=True)
    subset['max'] = subset.apply(lambda x: max(x.min(), x.max(), key=abs), axis=1)
    subset = subset.reindex(subset['max'].abs().sort_values(ascending=False, inplace=False).index)
    #give brief view of top and bottom vals
    print(subset['max'].head(n=6))
    print(subset['max'].tail(n=6))
    print(len(subset.index.values))
    if len(subset.index.values) > lim:
        result_list = subset.index.values[:lim]
    else:
        result_list = subset.index.values
    return result_list


def plot_ex_genes(df, h, n=10, mode='ind'):
    """ bundeled with find plot genes and plot concentrations to auto sample and plot dose-responsive
    gene plots from variety of larger magnitudes. will plot n (10) top mid and bottom genes from the list of those
    passing the threshold """
    gl = find_plot_genes(df)

    if len(gl) > 30:
        genes = []
        mid = round(len(gl)/2)
        genes.extend(gl[:n])
        genes.extend(gl[n-5:n+5])
        genes.extend(gl[-n:])
    else:
        genes = gl

    plot_concentrations(df, h, genes=genes, label=True, mode=mode)


def plot_concentrations(df, h, genes='test2', label=False, mode='ind', incr='dflt', outpath='dflt', fn='dflt', maxx='dflt',
                        test=False):
    """ plotting concentration plots on a per-gene basis from a df/header
    outpath for figures is optionally specified, genes can be passed in as
    a list, left as 'test' for a single gene, or be 'all'.
    the mode= ind,med,avg will either plot individual reps w/ same x value or combine
    reps together using either 'med' median or 'avg' average across reps

    assumes broken down by name and dose, and only within one batch
    """
    # parametetrs controlling the optional labels below each cohort
    txt_args = {'fontsize': 8,
                'rotation': 90,
                'fontweight': 'bold'}
    genes = gt.get_genes(genes, df=df)
    # define outpath directory, create if necessary
    if outpath is 'dflt':
        outpath = os.path.join(gt.check_desktop(), 'output_figs')
    try:
        os.mkdir(outpath)
    except:
        pass
    # define title
    if fn is not 'dflt':
        name = fn
    else:
        try:
            name = df.name
        except AttributeError:
            name = h.index[0].split(':')[0]
    # set the color pallet and spacing/sizing levels (figsize tuned to these)
    cmap = plt.get_cmap('tab10')
    if incr == 'dflt':
        incr = 10
    sincr = 20

    # sort the sample wells in desired order, by name and dose for test
    d, h = gt.dsub(df, h, {'type':'test'})
    df = d

    # create pert list for plot, strip batch if there's only one batch
    pert_list = []
    print(h['name'].unique())
    for n in h['name'].unique():
        pert_list.append('{}'.format(n))

    # if there are multiple reps adjust figure width to account
    # reps is for each name and dose combo, how many are there?
    #num_reps = round(h.groupby('name')['dose'].nunique().mean())
    ndoses = h.groupby('name')['dose'].nunique().max()
    nnames = h.name.nunique()
    print(name)
    print('num doses ', ndoses)
    print('name list ', len(pert_list))

    if isinstance(genes, str):
        genes = [genes]

    if maxx == 'dflt':
        # calc x range with length of vector corrected by reps, plus spacing btwn
        # basewidth = (len(d.iloc[0]) / num_reps) * incr
        # pert_buffer = (len(pert_list)) * 1 * incr
        pad = 8 * incr
        # maxx = basewidth + pert_buffer + pad
        maxx = (incr * nnames * ndoses) + (incr * 2 * nnames)

    for g in genes:
        # set initial color counters and x starting position
        ci = 0
        x_pos = 15
        # select vector for current gene
        dat = df.loc[g]
        # determine the max range of x axis
        maxv = round(max(abs(dat))) + 3
        ax = format_concentration_plot(maxx, maxy=maxv)
        ax.set_ylabel(g)
        mytitle = name + ' - ' + g
        print(mytitle)
        ax.set_title(mytitle)
        x_init = 0
        names = h['name'].apply(lambda x: str(x)).unique()
        for n in names:
            # increment through colors in cmap
            color = cmap(ci)
            ci += .1
            if ci > .9:
                ci = 0
            sub = h[h['name']==n]
            doses = sorted(sub['dose'].unique(), key=lambda x: float(x))
            sizes = [(x + 1) * sincr for x in range(len(doses))]
            for d, s in zip(doses, sizes):
                args = {'name': n, 'dose': d}
                wids = gt.hsub(h, args).index.values
                y_vals = dat[wids].values
                if mode == 'avg':
                    y_vals = np.mean(y_vals)
                if mode == 'med':
                    y_vals = np.median(y_vals)
                try:
                    x_vals = [x_pos] * len(y_vals)
                except TypeError:
                    x_vals = x_pos
                # plot the current vals with specified color and size
                ax.scatter(x_vals, y_vals, c=color, s=s)
                x_pos += incr
            # put spacing between perts
            if label is True:
                # n = ' '.join([n, d])
                x_label = (x_init + x_pos) / 2
                ax.text(x_label, -(maxv + 1), n, color=color, **txt_args)
            x_pos += (incr * 2)
            x_init = x_pos
        plt.savefig(os.path.join(outpath, mytitle + '.png'), bbox_inches='tight')
        plt.close()
        if test is True:
            print('test mode, exiting after one image')
            break


def plot_gene(sample_set, h=None, name='dflt', outpath='dflt', close=True, width=8):
    """ basic plot gene finction, if header is provided will apply color coding blue = veh, red = poscon """
    if name == 'dflt':
        name = sample_set.name
    if outpath == 'dflt':
        outpath = gt.dflt_outpath(fldr_name='dflt')
    xrange = len(list(sample_set))
    dtype = check_plottype(sample_set)
    #print('dtype is ', dtype)
    ax = format_concentration_plot(xrange, ptype=dtype, width=width)
    ax.scatter(range(xrange), sample_set.values, color='grey')
    ax.set_title(name)
    if h is not None:
        h['order'] = np.arange(1,len(h)+1)
        dv, hv = gt.dsub(sample_set, h, {'type': 'vehicle'})
        ax.scatter(hv.order, dv.values, color='blue')
        dp, hp = gt.dsub(sample_set, h, {'type': 'poscon'})
        ax.scatter(hp.order, dp.values, color='red')
    if close is True:
        plt.savefig(os.path.join(outpath, name + '.png'))
        plt.close()
    else:
        return ax


def make_dotplot(vctr, wdict=None, title='dflt', outpath='dflt', legend=False, width=5):
    """ passing in a series, label and well dictionary of highlighted cohorts with name: [wells] """
    if outpath == 'dflt':
        outpath = gt.dflt_outpath(fldr_name='output figs')
    cmap = plt.get_cmap('tab10')
    xrange = len(list(vctr))
    dtype = check_plottype(vctr.iloc[2])
    #print('dtype is ', dtype)
    ax = format_concentration_plot(xrange, ptype=dtype, width=width)
    # set additional title and axis
    ax.set_ylabel(vctr.name, fontsize=12)
    if title == 'dflt':
        title = vctr.index.values[0].split(':')[0] + ' - ' + vctr.name
    ax.set_title(title)
    awells = list(vctr.index)
    # plot primary data
    plt.plot(vctr.values, color='silver', marker='o', ls='', markersize=5, mew=0)
    ci = 0
    allwells = []
    [allwells.extend(w) for w in wdict.values()]
    if not any([w in vctr.index for w in allwells]):
        print('well dictionary not aligned, attempting patch')
        wdict2 = {}
        pname = vctr.index.values[0].split(':')[0]
        for name, wells in wdict.items():
            wdict2[name] = [pname + ':' + w for w in wells]
        wdict = wdict2
    if wdict is not None:
        mycolors, mynames = [], []
        for name, wells in wdict.items():
            color = cmap(ci)
            mycolors.append(color)
            mynames.append(name)
            ci += 1
            if ci > 9:
                ci = 0
            mywells = [x for x in wells if x in vctr.index]
            if len(mywells) == 0:
                print('no wells left')
            sety = vctr[mywells]
            setx = [awells.index(w) for w in mywells]
            plt.plot(setx, sety, color=color, marker='o', ls='', markersize=5, mew=0)
    if legend is not False:
        leg_dict = dict(zip(mynames, mycolors))
        # create a patch (proxy artist) for every color
        patches = [mpatches.Patch(color=mycolors[i], label=text) for i, text in enumerate(leg_dict.keys())]
        # put those patched as legend-handles into the legend
        lgd = plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    filename = title + '.png'
    plt.savefig(os.path.join(outpath, filename))
    plt.close()


def check_plottype(data):
    """ auto detect scale and processing level of data passed """
    try:
        if data.min() > 3:
            dtype = 'gct'
        elif data.mean() > 100:
            dtype = 'csv'
        else:
            dtype = 'zs'
    except ValueError:
        if any(data.min() > 3):
            dtype = 'gct'
        elif any(data.mean() > 100):
            dtype = 'csv'
        else:
            dtype = 'zs'
    #print(dtype)
    return dtype


def plot_cohorts(vdict, outpath='dflt', mode='sep', dtype='auto', maxx='dflt', title='foo', label=True, incr=1, size=20):
    """ given a dictionary with name : [values], plot them all on the same plot but with
    different color (and optional size). Can be promiscuity values, or concentration range
    breakdown will give a dictionary in return, optionally
    dtype 'auto' will figure out either gct or zs, but can be specified """
    if outpath == 'dflt':
        outpath = gt.dflt_outpath()
    # set the color pallet and spacing/sizing levels (figsize tuned to these)
    cmap = plt.get_cmap('tab10')
    # parametetrs controlling the optional labels below each cohort
    txt_args = {'fontsize': 8,
                'rotation': 90,
                'fontweight': 'bold'}
    # set initial color counters and x starting position
    ci = 0
    x_pos = 1
    # is this duplicated with the stuff below?
    try:
        maxv = round(max([max(abs(max(v)), abs(min(v))) for v in vdict.values()]))
    except:
        try:
            maxv = round(max([max(abs(v)) for v in vdict.values()]))
        except:
            try:
                maxv = round(max([abs(v).max() for v in vdict.values()]))
            except:
                maxv = round(max([abs(v).max().max() for v in vdict.values()]))
    maxv += maxv*0.1
    # calc x range with length of vector corrected by reps, plus spacing btwn
    if mode is 'sep':
        try:
            maxx = sum([len(x.columns) for x in vdict.values()]) * incr + incr
        except:
            maxx = sum([len(x) for x in vdict.values()]) * incr + incr
    elif mode is 'tog':
        maxx = len(vdict.keys()) * incr + incr
    # pull out the first value set to check for plot formatting
    for i, vals in enumerate(vdict.values()):
        if i >= 1:
            break
        myvals = vals
    # create and baseline format plot
    if dtype == 'auto':
        dtype = check_plottype(myvals)
    # adjust plot type if auto adjusted
    if dtype == 'zs' and maxv > 10:
        maxy = round(maxv + 1)
        ax = format_concentration_plot(maxx, ptype=dtype, maxy=maxy)
    else:
        ax = format_concentration_plot(maxx, ptype=dtype)
    ax.set_xlim([0,maxx + 1])
    ax.set_xlabel('')
    # determine title,
    if title is 'foo':
        try:
            title = myvals.name
        except:
            try:
                title = myvals.columns[0]
            except:
                title = 'foo'
    ax.set_title(title)

    if dtype == 'gct':
        y_label = min(ax.get_ylim()) * 0.75
    elif dtype == 'zs':
        y_label = min(ax.get_ylim()) * 1.2

    for n, vals in vdict.items():
        #try:
        #    vals = vals.values[0]
        #except:
        #    print('vals error...')
        #    pass
        # increment through colors in cmap
        color = cmap(ci)
        ci += 1
        if ci > 9:
            ci = 0
        # catch to handle empty values
        if len(vals) == 0:
            xlength = 1
            vals = []
        else:
            xlength = len(vals)
        # set x coordinates for values
        if mode is 'sep':
            x_vals = [x_pos + (x*incr) for x in range(xlength)]
            x_pos = max(x_vals) + incr
        elif mode is 'tog':
            x_vals = [x_pos] * xlength
            x_pos += incr
        # plot the current vals with specified color and size
        # implement catches in case no values for a given entry
        try:
            ax.scatter(x_vals, vals, color=color, s=size)
        except:
            ax.scatter(x_vals, [0], color='white', s=size)
        #print(vals)
        # then add label for each cohort below
        if label is True:
            if len(x_vals) > 1:
                x_label = (x_vals[0] + x_vals[-1])/ 2
            else:
                x_label = x_vals[0] - 1
            ax.text(x_label, y_label, n, color=color, **txt_args)
    #return ax
    plt.savefig(os.path.join(outpath, title + '.png'), bbox_inches='tight')
    plt.close()


def make_plateplot(vctr, name='dflt', outpath='dflt', label='dflt', cmap='inferno'):
    """ plot feature in 384-well plate format """
    vctr = prep_vctr(vctr)
    plot_plateplot(vctr, name=name, outpath=outpath, label=label, cmap=cmap)


def prep_vctr(vctr, order='row'):
    """ format and ensure vector is correct length and order """
    new_index = gt.get_awells()
    if any(vctr.index.str.contains(':')):
        pname = vctr.index.values[0].split(':')[0]
        new_index = [pname + ':' + x for x in new_index]
    awells = pd.Series(np.nan * 384, index=new_index)
    # merge the new wells onto an array of 384 in order to keep spacing
    fvctr = awells.combine_first(vctr)
    #fvctr = pd.concat([awells,vctr], ignore_index=True)
    if order is 'col':
        print('plotting in column order')
        # create a new sorted vector by sorting by number first, then row
        svctr = fvctr.loc[sorted(fvctr.index.values, key=lambda x: (x[1:],x[0]))]
    else:
        # sort vector samples by well row value by default
        svctr = fvctr.loc[sorted(fvctr.index)]
    svctr.name = vctr.name
    return svctr


def plot_plateplot(vctr, name='dflt', outpath='dflt', label='dflt', cmap='inferno', ncats=None, clrbar=True):
    """ will plot a 384 well plate with the values passed in the Series object vector, will 
    map to all wells and plot with auto-adjusted colors in the provided map with colorbar w/
    values if clrbar is True. Otherwise can pass dictionary into the clrbar variable to 
    plot a separate individual legend with keys as the name and values as the converted integer used
    to plot the map """
    if name == 'dflt':
        name = vctr.index.values[0].split(':')[0] + '-' + vctr.name
    # elif '_' not in name:
    #     name = name + ' - ' + vctr.name
    else:
        name = name
    if outpath == 'dflt':
        outpath = gt.dflt_outpath()
    if label == 'dflt':
        #label = vctr.name
        label = name
    fig, ax = plt.subplots()
    # set additional title and axis
    ax.set_title(label, y=1.1, fontsize=16)
    row_labels = list(string.ascii_uppercase[0:16])
    ax.set_yticks(list(np.arange(16)))
    ax.set_yticklabels(row_labels, fontsize=8)
    col_labels = list(np.arange(1,25))
    ax.set_xticks(list(np.arange(0,24)))
    ax.set_xticklabels(col_labels, fontsize=9)
    ax.tick_params(labelright=True, labeltop=True)
    # this sets the tick length to zero, but leaves labels
    plt.tick_params(axis=u'both', which=u'both',length=0)
    # reshape array and plot
    try:
        d = vctr.values.reshape(16,24)
    except:
        print('error in reshape')
        return
    if ncats is not None:
        im = plt.imshow(d, interpolation='nearest', cmap=cmap, vmin=0, vmax=ncats)
    else:
        im = plt.imshow(d, interpolation='nearest', cmap=cmap)
    # use matplotlib axes1 to keep colorbars in line with figs
    if clrbar is True:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.3)
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.tick_params(labelsize=9)
        # simplify colorbar to 5 points including max/min
        mx, mn = vctr.max(), vctr.min()
        mid = (mx + mn) / 2
        svth = mid + ((mx - mid)/2)
        twth = (mid - ((mx - mid)/2))
        things = [mn, twth, mid, svth, mx]
        # things = [mn, mid, mx]
        thingsl = ['{:.1f}'.format(x) for x in things]
        cbar.set_ticks(things)
        cbar.set_ticklabels(thingsl)
    elif clrbar is not True:
        # get the colors of the values, according to the 
        # colormap used by imshow
        leg_dict = clrbar
        colors = [im.cmap(im.norm(value)) for value in leg_dict.values()]
        # create a patch (proxy artist) for every color 
        patches = [ mpatches.Patch(color=colors[i], label=text ) for i, text in enumerate(leg_dict.keys()) ]
        # put those patched as legend-handles into the legend
        lgd = plt.legend(handles=patches, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0. )
    #plt.tight_layout()
    #if not outpath.endswith('.png'):
    #    outpath += '.png'
    outpath = os.path.join(outpath, name + '.png')
    try:
        fig.savefig(outpath, bbox_extra_artists=(lgd,), bbox_inches='tight')
    except:
        plt.savefig(outpath)
    plt.close()



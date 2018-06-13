#!/usr/bin/python


import gct, os, pt, gt, string, scipy, math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import collections as cll
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib import rcParams


def plot_landmark_concs(df, h, maxy=12, genes='test100', outpath='dflt', test=False):
    """ plot many or all landmarks, should pass in a subset dataframe and header which
    should be the consensus ZS file. can contain many different names + doses, will auto breakdown by 'nd'
    a single line per gene is plotted for the ZS across all concentrations """
    # txt_args = {'fontsize': 8,
    #             'rotation': 90,
    #             'fontweight': 'bold'}

    if outpath is 'dflt':
        outpath = gt.dflt_outpath(fn=gt.get_shn(df.name))
    for ds, hs in gt.breakdown(df, h, 'n', dic=False):
        xrange = len(hs.dose.unique())-2
        ax = format_concentration_plot(xrange, maxy=maxy, width=4)
        ax.tick_params(axis='x', bottom='on', top='off', labelbottom='on')
        ax.set_xticks(range(len(hs.dose.unique())))
        nname = hs['name'][0]
        ax.set_title(ds.name[:6] + ' - ' + nname, fontsize=14)
        for g in gt.get_genes(genes, df=df):
            data = ds.loc[g,:]
            ax.plot(data.values, linewidth=0.3)
        # if label is True:
        #     x_label = (x_vals[0] + x_vals[-1])/ 2
        #     ax.text(x_label, y_label, n, color=color, **txt_args)
        plt.savefig(outpath + '-' + nname + '.png')
        plt.close()
        if test is True:
            print('stopping after one iteration')
            break


def get_euclidean(df):
    Y = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(df.T, 'euclidean'))
    return Y


def plot_euclidean(df, labels, rot=None, fontsize=None, tick_denom=1, dat=None):
    ''' takes dataframe and plots square euclidean distance, labels w/ labels
    must define ax = plot_euclidean() to get it back an operable object'''
    Y = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(df.T, 'euclidean'))
    #if up is not None:
    #    Y = Y.clip(max=up)
    if dat is not None:
        Y = dat
    ax = sns.heatmap(Y, square=True, cmap='bwr_r')
    # set automatic tick positioning
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tick_denom))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tick_denom))
    # insert blanks start + end of labels to get it to display right, a little odd but works
    labels = list(labels)
    labels.insert(0,'')
    labels.insert(0,'')
    labels.append([''])
    if rot is None:
        rot = 90
        hor_align = 'center'
    else:
        hor_align = 'right'
    if fontsize == None:
        #ax.set_xticklabels(labels, rotation=40, ha='right')
        ax.set_xticklabels(labels, rotation=rot, ha=hor_align)
        #ax.set_yticklabels(labels[::-1], rotation=0, va='center')
        ax.set_yticklabels(labels, rotation=0, va='center')
    else:
        #ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=fontsize)
        ax.set_xticklabels(labels, rotation=rot, ha=hor_align, fontsize=fontsize)
        #ax.set_yticklabels(labels[::-1], rotation=0, va='center', fontsize=fontsize)
        ax.set_yticklabels(labels, rotation=0, va='center', fontsize=fontsize)
    plt.tight_layout()
    return ax


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



## under construction !!
def format_barview_plot(enrich, hlight_idx):
    """ plot a barview, with edf passed in or just  """
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


def find_plot_genes(df, thresh=8, lim=50):
    """ open a gct file and make list of genes which have ZS above 8 to plot concentrations """
    print(df.name)
    subset = df[abs(df) > thresh]
    subset.dropna(axis=0, how='all', inplace=True)
    subset['max'] = subset.apply(lambda x: max(x.min(), x.max(), key=abs), axis=1)
    subset.reindex(subset['max'].abs().sort_values(ascending=False, inplace=False).index)
    print(subset['max'].head(n=10))
    print(len(subset.index.values))
    if len(subset.index.values) > lim:
        result_list = subset.index.values[:lim]
    else:
        result_list = subset.index.values
    return result_list


def plot_concentrations(f, genes='test2', label=False, mode='ind', incr='dflt', outpath='dflt', fn='dflt', maxx='dflt', test=False):
    """ plotting concentration plots on a per-gene basis from a gct file
    outpath for figures is optionally specified, genes can be passed in as
    a list, left as 'test' for a single gene, or be 'all'.
    the mode= ind,med,avg will either plot individual reps w/ same x value or combine
    reps together using either 'med' median or 'avg' average across reps
    -- now also accepts [d,h] in place of filepath as  primary argument
    assumes broken down by name and dose... """
    if isinstance(f, str):
        df, h = gct.extractgct(f)
        print(df.name)
    else:
        try:
            df, h = f[0], f[1]
        except (IndexError, KeyError):
            exit(code='problem with plot concentration input, needs file or [df,h]')


    # parametetrs controlling the optional labels below each cohort
    txt_args = {'fontsize': 8,
                'rotation': 90,
                'fontweight': 'bold'}

    genes = gt.get_genes(genes, df=df)

    # define outpath directory, create if necessary
    if outpath is 'dflt':
        outpath = os.path.join(pt.check_desktop(), 'output_figs')
    try:
        os.mkdir(outpath)
    except:
        pass

    # set the color pallet and spacing/sizing levels (figsize tuned to these)
    cmap = plt.get_cmap('tab10')

    if incr == 'dflt':
        incr = 10
    sincr = 20

    # sort the sample wells in desired order
    d, h = gt.dsub(df, h, {'type':'test'})

    ord_wells = h.sort_values(['name', 'dose']).index.values
    df = d[ord_wells]

    if fn is not 'dflt':
        name = fn
    else:
        try:
            name = df.name
        except AttributeError:
            name = h.addr[0].split(':')[0]

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
        names = sorted(h['name'].apply(lambda x: str(x)).unique())
        for n in names:
            # increment through colors in cmap
            color = cmap(ci)
            ci += 1
            if ci > 9:
                ci = 0
            sub = h[h['name']==n]
            doses = sorted(sub['dose'].unique(), key=lambda x: float(x))
            sizes = [(x + 1) * sincr for x in range(len(doses))]
            for d, s in zip(doses, sizes):
                args = {'name':n, 'dose':str(d)}
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


def plot_gene(sample_set):
    xrange = len(list(sample_set))
    dtype = check_plottype(sample_set.iloc[0])
    print('dtype is ', dtype)
    ax = format_concentration_plot(xrange, ptype=dtype)
    ax.scatter(range(xrange), sample_set.values)
    ax.set_title(sample_set.name)
    outpath = gt.dflt_outpath('dflt', 'output_figs')
    plt.savefig(os.path.join(outpath, sample_set.name + '.png'))


def plot_gene_wtypes(vctr, h):
    xrange = len(list(vctr.values))
    dtype = check_plottype(vctr.iloc[0])
    print('dtype is ', dtype)
    ax = format_concentration_plot(xrange, ptype=dtype)
    ax.scatter(range(xrange), vctr.values, color='silver')

    ax.set_title(vctr.name)
    outpath = gt.dflt_outpath('dflt', 'output_figs')
    plt.savefig(os.path.join(outpath, vctr.name + '.png'))


def make_dotplot(vctr, wdict=None, title='dflt', outpath='dflt', legend=False):
    """ passing in a series, label and well dictionary of highlighted cohorts with name: [wells] """
    if outpath == 'dflt':
        outpath = gt.dflt_outpath(fldr_name='output figs')
    cmap = plt.get_cmap('tab10')
    xrange = len(list(vctr))
    dtype = check_plottype(vctr.iloc[2])
    #print('dtype is ', dtype)
    ax = format_concentration_plot(xrange, ptype=dtype)
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
    print(dtype)
    return dtype


def plot_cohorts(vdict, outpath='dflt', mode='sep', maxx='dflt', title='foo', label=True, incr=1, size=20):
    """ given a dictionary with name : [values], plot them all on the same plot but with
    different color (and optional size). Can be promiscuity values, or concentration range
    breakdown will give a dictionary in return, optionally"""
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
    # determine the max range of x axis
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
            maxx = sum([len(x.columns) for x in vdict.values()]) * incr
        except:
            maxx = sum([len(x) for x in vdict.values()]) * incr
    elif mode is 'tog':
        maxx = len(vdict.keys()) * incr
    print('max x is ', maxx)

    # pull out the first value set to check for plot formatting
    for i, vals in enumerate(vdict.values()):
        if i >= 1:
            break
        myvals = vals
    # determine title
    if title is 'foo':
        try:
            title = myvals.name
        except:
            try:
                title = myvals.columns[0]
            except:
                title = 'foo'
    # create and baseline format plot
    dtype = check_plottype(myvals)
    ax = format_concentration_plot(maxx, ptype=dtype)
    ax.set_xlim([0,maxx + 1])
    ax.set_xlabel('')
    ax.set_title(title)

    y_label = min(ax.get_ylim()) * 0.75

    for n, vals in vdict.items():
        try:
            vals = vals.values[0]
        except:
            print('vals error...')
            pass
        # increment through colors in cmap
        color = cmap(ci)
        ci += 1
        if ci > 9:
            ci = 0
        if mode is 'sep':
            x_vals = [x_pos + (x*incr) for x in range(len(vals))]
            x_pos = max(x_vals) + incr
        elif mode is 'tog':
            x_vals = [x_pos] * len(vals)
            x_pos += incr
        # plot the current vals with specified color and size
        ax.scatter(x_vals, vals, color=color, s=size)
        print(vals)
        # then add label for each cohort below
        if label is True:
            x_label = (x_vals[0] + x_vals[-1])/ 2
            ax.text(x_label, y_label, n, color=color, **txt_args)
    #return ax
    plt.savefig(os.path.join(outpath, title + '.png'), bbox_inches='tight')
    plt.close()


# incomplete!!
def plot_promiscuity(path):
    d = pd.read_table(path, index_col=0)
    sort_cols = ['name', 'dose']
    if 'cell' in d.columns:
        sort_cols.insert(0, 'cell')
    d.sort_values(sort_cols, inplace=True)
    sort_dict = cll.OrderedDict()
    plot_dict = cll.OrderedDict()
    for c in sort_cols:
        sort_dict[c] = d[c].unique()


def make_plateplot(vctr, name='dflt', outpath='dflt', label='dflt', cmap='inferno'):
    vctr = prep_vctr(vctr)
    plot_plateplot(vctr, name=name, outpath=outpath, label=label, cmap=cmap)


def prep_vctr(vctr, order='row'):
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
        name = vctr.name
    else:
        name = name + ' - ' + vctr.name
    if outpath == 'dflt':
        outpath = gt.dflt_outpath(fn=name)
    if label == 'dflt':
        label = vctr.name
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
    if not outpath.endswith('.png'):
        outpath += '.png'
    try:
        fig.savefig(outpath, bbox_extra_artists=(lgd,), bbox_inches='tight')
    except:
        plt.savefig(outpath)
    plt.close()


def main():
    plot_concentrations('/Users/WRB/Desktop/data/DEV098_ZSVCQNORM_n330x978.gct', genes='test1')
    # pass

if __name__ == '__main__':
    main()

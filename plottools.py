#!/usr/bin/python

import gct, os, pt, gt
import matplotlib.pyplot as plt
import numpy as np


def plot_euclidean(df, labels, fontsize=None):
  ''' takes dataframe and plots square euclidean distance, labels w/ labels
  must define ax = plot_euclidean() to get it back an operable object'''
  Y = scipy.spatial.distance.squareform(scipy.spatial.distance.pdist(df.T, 'euclidean'))
  ax = sns.heatmap(Y, square=True, cmap='bwr_r')
  if fontsize == None:
    ax.set_xticklabels(labels, rotation=40, ha='right')     
    ax.set_yticklabels(labels[::-1], rotation=0, va='center')
  else:
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=fontsize)     
    ax.set_yticklabels(labels[::-1], rotation=0, va='center', fontsize=fontsize)
  plt.tight_layout()
  return ax


def fix_ticklabels(ax, l):
    ax.set_xticklabels(l, rotation=40, ha='right')
    #ax.set_yticklabels(l[::-1], rotation=0, va='center')
    ax.set_yticklabels(l, rotation=0, va='center')
    return ax


def set_titles(ax, main, xlabel, ylabel, xrot=0, yrot=0):
    """ shortcut to label axis with fig title, axis labels and set x+y tickmark rotation """
    if main == None:
        sys.exit('titles(ax, main, xtitle, ytitle, rot=0)')

    ax.set_title(main)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.set_xticklabels(ax.xaxis.get_majorticklabels(), rotation=xrot)
    ax.set_yticklabels(ax.yaxis.get_majorticklabels(), rotation=yrot)
    

def format_concentration_plot(maxx, ptype='zs', maxy=10):
    """ basic simple y-axis only plot formatting. set the ylims externally
    (so this is flexible for zs and gct) """
    width = maxx / 10
    if width > 12:
        width = 12
    fig = plt.figure(figsize=[width, 4])
    ax = fig.add_subplot(111)
    ax.set_facecolor("white")
    # format graph
    if ptype is 'zs':
        ax.set_ylim([-maxy, maxy])
        ax.set_yticks([-maxy, 0, maxy])
    elif ptype is 'gct':
        ax.set_ylim([2, 16])
        ax.set_yticks(np.arange(2, 18, 2))
    elif ptype is 'csv':
        ax.set_ylim([2, 16])
        if maxy > 500:
            ax.set_yticks(np.arange(0, 10000, 2))
        else:
            ax.set_yticks(np.arange(0, 10000, 2))
    ax.set_xlim([0, maxx+20])

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


def plot_concentrations(f, outpath='dflt', genes='test2', mode='ind'):
    """ plotting concentration plots on a per-gene basis from a gct file
    outpath for figures is optionally specified, genes can be passed in as
    a list, left as 'test' for a single gene, or be 'all'.
    the 'mode' will either plot individual reps w/ same x value or combine
    reps together using either 'med' median or 'avg' average across """

    df, h = gct.extractgct(f)
    print(df.name)

    # check for gene arguments, and fill 'test' and 'all'
    if genes == 'test1':
        genes = ['200678_x_at']
    elif genes == 'test2':
        genes = ['121_at', '219888_at', '218245_at', '206501_x_at', '203154_s_at',
                 '201614_s_at', '209682_at', '202324_s_at', '209603_at',
                 '200060_s_at', '202123_s_at', '201579_at']
    elif genes == 'all':
        genes = df.index.values

    # define outpath directory, create if necessary
    if outpath == 'dflt':
        outpath = os.path.join(pt.check_desktop(), 'output_figs')
    try:
        os.mkdir(outpath)
    except:
        pass

    # set the color pallet and spacing/sizing levels (figsize tuned to these)
    cmap = plt.get_cmap('tab10')
    incr = 10
    sincr = 20

    # sort the sample wells in desired order
    ord_wells = h[h['type'] == 'test'].sort_values(['batch', 'name', 'dose']).index.values
    dord = df[ord_wells]
    name = df.name

    names = sorted(h[h['type'] == 'test']['name'].unique())
    batches = h['batch'].dropna().unique()

    # create pert list for plot, strip batch if there's only one batch
    pert_list = []
    for b in batches:
        for n in names:
            pert_list.append('{}-{}'.format(n,b))
    if len(batches) == 1:
        pert_list = [x.rstrip('-A') for x in pert_list]
    print(pert_list)

    # if there are multiple reps adjust figure width to account
    dummy = pt.hsub(h, {'batch':'A', 'name':names[0]})
    dummy_doses = list(dummy['dose'].unique())
    num_reps = len(dummy[dummy['dose']==dummy_doses[0]].index.values)

    for g in genes:
        # set initial color counters and x starting position
        ci = 0
        x_pos = 15
        # select vector for current gene
        dat = dord.loc[g]
        # determine the max range of x axis
        maxv = round(max(abs(dat))) + 1
        # calc x range with length of vector corrected by reps, plus spacing btwn
        maxx = ((len(dat) * incr)/num_reps) + (len(pert_list)*2*incr) + 30
        ax = format_concentration_plot(maxv, maxx)
        ax.set_ylabel(g)
        title = name + ' - ' + g
        ax.set_title(title)
        for b in batches:
            for n in names:
                # increment through colors in cmap
                color = cmap(ci)
                ci += 1
                if ci > 9:
                    ci = 0
                sub = h[h['name']==n]
                doses = sorted(sub['dose'].unique())
                sizes = [(x + 1) * sincr for x in range(len(doses))]
                for d, s in zip(doses, sizes):
                    args = {'batch':b, 'name':n, 'dose':d}
                    wids = pt.hsub(h, args).index.values
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
                x_pos += (incr*2)
        plt.tight_layout()
        plt.close()
        plt.savefig(os.path.join(outpath, title + '.png'))


def plot_gene(sample_set):
    xrange = len(list(sample_set))
    if sample_set.iloc[0].min() > 3:
        dtype = 'gct'
    elif sample_set.iloc[0].mean() > 100:
        dtype = 'csv'
    else:
        dtype = 'zs'
    print(dtype)
    ax = format_concentration_plot(xrange, ptype=dtype)
    ax.scatter(range(xrange), sample_set.values)
    ax.set_title(sample_set.name)
    outpath = gt.dflt_outpath('dflt', 'output_figs')
    plt.savefig(os.path.join(outpath, sample_set.name + '.png'))


def plot_gene_cohorts(f, outpath='dflt', genes='test2', mode='ind'):
    """ plotting concentration plots on a per-gene basis from a gct file
    outpath for figures is optionally specified, genes can be passed in as
    a list, left as 'test' for a single gene, or be 'all'.
    the 'mode' will either plot individual reps w/ same x value or combine
    reps together using either 'med' median or 'avg' average across """

    df, h = gct.dfsubset(f, wids)
    print(df.name)

    # check for gene arguments, and fill 'test' and 'all'
    if genes == 'test1':
        genes = ['200678_x_at']
    elif genes == 'test2':
        genes = ['121_at', '219888_at', '218245_at', '206501_x_at', '203154_s_at',
                 '201614_s_at', '209682_at', '202324_s_at', '209603_at',
                 '200060_s_at', '202123_s_at', '201579_at']
    elif genes == 'all':
        genes = df.index.values

    # define outpath directory, create if necessary
    if outpath == 'dflt':
        outpath = os.path.join(pt.check_desktop(), 'output_figs')
    try:
        os.mkdir(outpath)
    except:
        pass

    # set the color pallet and spacing/sizing levels (figsize tuned to these)
    cmap = plt.get_cmap('tab10')
    incr = 10
    sincr = 20

    # sort the sample wells in desired order
    ord_wells = h[h['type'] == 'test'].sort_values(['batch', 'name', 'dose']).index.values
    dord = df[ord_wells]
    name = df.name

    names = sorted(h[h['type'] == 'test']['name'].unique())
    batches = h['batch'].dropna().unique()

    # create pert list for plot, strip batch if there's only one batch
    pert_list = []
    for b in batches:
        for n in names:
            pert_list.append('{}-{}'.format(n,b))
    if len(batches) == 1:
        pert_list = [x.rstrip('-A') for x in pert_list]
    print(pert_list)

    # if there are multiple reps adjust figure width to account
    dummy = pt.hsub(h, {'batch':'A', 'name':names[0]})
    dummy_doses = list(dummy['dose'].unique())
    num_reps = len(dummy[dummy['dose']==dummy_doses[0]].index.values)

    for g in genes:
        # set initial color counters and x starting position
        ci = 0
        x_pos = 15
        # select vector for current gene
        dat = dord.loc[g]
        # determine the max range of x axis
        maxv = round(max(abs(dat))) + 1
        # calc x range with length of vector corrected by reps, plus spacing btwn
        maxx = ((len(dat) * incr)/num_reps) + (len(pert_list)*2*incr) + 30
        ax = format_concentration_plot(maxv, maxx)
        ax.set_ylabel(g)
        title = name + ' - ' + g
        ax.set_title(title)
        for b in batches:
            for n in names:
                # increment through colors in cmap
                color = cmap(ci)
                ci += 1
                if ci > 9:
                    ci = 0
                sub = h[h['name']==n]
                doses = sorted(sub['dose'].unique())
                sizes = [(x + 1) * sincr for x in range(len(doses))]
                for d, s in zip(doses, sizes):
                    args = {'batch':b, 'name':n, 'dose':d}
                    wids = pt.hsub(h, args).index.values
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
                x_pos += (incr*2)
        plt.tight_layout()
        plt.savefig(os.path.join(outpath, title + '.png'))
        plt.close()


def prep_vctr(vctr, order='row'):
    new_index = get_awells()
    if any(vctr.index.str.contains(':')):
        pname = vctr.index.values[0].split(':')[0]
        new_index = [pname + ':' + x for x in new_index]
    awells = pd.Series(np.nan * 384, index=new_index)
    # merge the new wells onto an array of 384 in order to keep spacing
    fvctr = awells.combine_first(vctr)
    if order is 'col':
        print('plotting in column order')
        # create a new sorted vector by sorting by number first, then row
        svctr = fvctr.loc[sorted(fvctr.index.values, key=lambda x: (x[1:],x[0]))]
    else:
        # sort vector samples by well row value by default
        svctr = fvctr.loc[sorted(fvctr.index)]
    return svctr


def make_plateplot(vctr, outpath='dflt', label=vctr.name, cmap='inferno'):
    if outpath is 'dflt':
        outpath = gt.dflt_outpath(fn=vctr.name)
    fig, ax = plt.subplots()
    # set additional title and axis
    ax.set_ylabel(label, fontsize=16)
    row_labels = list(string.ascii_uppercase[0:16])
    row_range = list(np.arange(16))
    ax.set_yticks(row_range, row_labels, fontsize=8)
    col_labels = list(np.arange(1,25))
    col_range = list(np.arange(0,24))
    ax.set_xticks(col_range, col_labels, fontsize=9)
    ax.set_tick_params(labelright=True, labeltop=True)
    # this sets the tick length to zero, but leaves labels
    plt.tick_params(axis=u'both', which=u'both',length=0)
    # reshape array and plot
    d = vctr.values.reshape(16,24)
    im = plt.imshow(d, interpolation='nearest', cmap=cmap)
    # use matplotlib axes1 to keep colorbars in line with figs
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
    plt.tight_layout()
    if not outpath.endswith('.png'):
        outpath += '.png'
    plt.savefig(outpath)


def main():
    pass

if __name__ == '__main__':
    main()

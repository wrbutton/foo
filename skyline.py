#!/usr/bin/python

'''
1/28/17 - separeted out from the growing plateanalysis script
into justh the skyline plotting script, calling functions
from plateanalysis as needed
'''
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import gct, os, gt, pa


def make_plots_from_list(d, h, welllist, cats='nd', outpath='dflt', test=False):
    """ take a list of well ids (can be 3 char) and gets the matching reps with provided
    cats, and then plots skylines from the reps """
    pname = d.name
    if outpath == 'dflt':
        outpath = gt.dflt_outpath()
    for w in welllist:
        wells = gt.get_well_reps(h, w, cats, df=True)
        try:
            if wells == 'empty':
                continue
        except ValueError:
            pass
        name = wells.iloc[0]['name'] + '-' + wells.iloc[0]['dose']
        name = name.replace('.', ',')
        title = pname + '-' + name
        wids = wells.index.values
        wids = [x for x in wids if x in d.columns]
        myoutpath = os.path.join(outpath, title)
        new_skyline(d[wids], title=title, outpath=myoutpath)
        if test is True:
            print('test mode, exiting after one image')
            break

def sep_updown(well, t=0):
    """ pass a zscore instance and get back separate positive and negative vectors
    the t threshold argument nulls out values below the indicated threshold """
    wup = well.copy()
    wup[wup < t] = np.nan
    wdn = well.copy()
    wdn[wdn > -t] = np.nan
    return wup, wdn


def well_consensus(wells, name):
    """ takes list of pandas Series of Z-scored data as inputs, will return consensus
    of the minimum absolute score per gene """
    consensus = []
    wl = len(wells)
    for g in wells[0].index:
        vals = [wells[i][g] for i in range(wl)]
        if all([v > 0 for v in vals]):
            consensus.append(round(min(vals),3))
        elif all([v < 0 for v in vals]):
            consensus.append(round(max(vals),3))
        else:
            consensus.append(0)
    consensus = pd.Series(consensus, index=wells[0].index, name=name)
    return consensus


def format_skyline_plot(ax):
    ax.set_facecolor("white")
    # turn off the axes boxes except for y axis
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(True)
    # set the subplot title
    ax.yaxis.set_label_position('right')
    ax.tick_params(axis='x', bottom='off', top='off', labelbottom='off')
    ax.tick_params(axis='y', right='off', left='on')
    ax.tick_params(axis='y', direction='in', length=3, width=1)
    return ax


# add support for dictionary of highlights to have different classes/colors of highlight genes
# with a labeled legend
def plot_skyline(wells, highlights=None, labels=None, cons=True, outpath='dflt', t=1, maxv=10, title='auto', line1=None, line2=None):
    """ well can either be a single instance (pd Series) or a list of replicate instances, or dataframe,
     in which case the consensus skyline will be generated and plotted. maxv is yscale, 10 dflt or 'auto'
     which sets based upon max value in the dataset. line 1 and line 2 may contain metadata
     for labelling the plot, and highlights accepts a list of genes to be accented
     cons (consensus) default is True, and will assemble local consensus from passed wells.
     cons can also be set to False or None, where just the individual wells will be plotted.
     lastly, a precomputed consensus vector can be passed in which will be used instead of calculating.
     labels will add a second line of text to the right side y label of each well """
    # check if wells is long its just a single series, and put it in a list container
    # otherwise if wells is multiple replicate instances (length 2-5) leave as is
    try:
        if len(wells.columns) > 1:
            wells = [wells.loc[:,i] for i in wells.columns]
    except AttributeError:
        if len(wells) > 10:
            wells = [wells]
        if len(wells) == 0:
            print('empty wells')
    if outpath is 'dflt':
        outpath = gt.dflt_outpath(fldr_name='output figs')
        if title is 'auto':
            fn = wells[0].name.replace(':','-')
        else:
            fn = title
        outpath = os.path.join(outpath, fn)
    # sort genes in ascending alpha order for consistency
    [w.sort_index(inplace=True) for w in wells]
    names, data = [], []
    for i,w in enumerate(wells):
        data.append(sep_updown(w, t=t))
        if labels is not None:
            names.append(w.name + '\n' + labels[i])
        else:
            names.append(w.name)
    # pass to external function to get consensus of wells if multiple instances
    if len(wells) > 1:
        if cons is True:
            c = well_consensus(wells, 'consensus')
            data.append(sep_updown(c, t=t))
            names.append('consensus')
        elif cons is None:
            pass
        elif cons is False:
            pass
        else:
            c = cons
            data.append(sep_updown(c, t=t))
            names.append('consensus')
    # define plot y axis max/min zscore by max value in data
    if maxv=='auto':
        maxval = 0
        for d in data:
            val1 = max(abs(d[0]))
            val2 = max(abs(d[1]))
            submax = max(val1, val2)
            if submax > maxval:
                maxval = submax
        maxv = round(maxval) + 2
    # create figure
    l = len(data)
    fig = plt.figure(figsize=(10, (l*2)+2), facecolor="white")
    for i, (upv, dnv) in enumerate(data):
        # define subplot in grid, counting down from the max length
        ax = fig.add_subplot(l, 1, l-i, facecolor="white")
        ax.axis([0, 978, -maxv, maxv])
        ax = format_skyline_plot(ax)
        ax.set_ylabel(names[i])
        ax.set_yticks([-maxv, 0, maxv])
        # plot up + down values in green/red
        ax.bar(range(978), upv.values, width=2, facecolor='green', edgecolor='green')
        ax.bar(range(978), dnv.values, width=2, facecolor='firebrick', edgecolor='firebrick')
        # plot optional highlight genes
        if highlights:
            # check if highlight is a string, if so tuck in list so looping works
            if isinstance(highlights, str):
                highlights = [highlights]
            for g in highlights:
                try:
                    vx = np.where(wells[0].index.values==g)[0][0]
                except:
                    vx = np.where(wells[0].index.values==g)[0]
                ax.bar(vx, maxv, width=2, facecolor='lightgray', edgecolor='lightgrey')
                ax.bar(vx, -maxv, width=2, facecolor='lightgray', edgecolor='lightgrey')
                ax.bar(vx, upv[g], width=2, facecolor='lime', edgecolor='lime')
                ax.bar(vx, dnv[g], width=2, facecolor='red', edgecolor='red')
        # white line at zero
        ax.bar(range(978), np.zeros(978), width=1, facecolor='grey', edgecolor='grey')

    # load in metadata from optional line arguments and incorporate into title on multiple lines
    if title is 'auto':
        title = pt.get_shn(outpath)
    if line1 is not None:
        if line2 is not None:
            title = title + '\n' + str(line1) + '\n' + str(line2)
        else:
            title = title + '\n' + str(line1)

    plt.suptitle(title)
    plt.savefig(outpath, facecolor="white")
    plt.close()



def new_skyline(wells, highlights=None, outpath='dflt', t=1, maxv=10, title='auto', info=None):
    # updated with pandas operators for easier consensus derivation
    """ well can either be a single instance (pd Series) or a list of replicate instances, or dataframe,
     in which case the consensus skyline will be generated and plotted. maxv is yscale, 10 dflt or 'auto'
     which sets based upon max value in the dataset. line 1 and line 2 may contain metadata
     for labelling the plot, and highlights accepts a list of genes to be accented """
    # assumes dataframe data type, or series
    if isinstance(outpath, str) and outpath is 'dflt':
        outpath = gt.dflt_outpath(fldr_name='output figs')
        if title == 'auto':
            try:
                fn = wells.iloc[:,0].name.replace(':','-')
            except:
                fn = wells.name.replace(':','-')
        else:
            fn = title
        outpath = os.path.join(outpath, fn)
    # sort genes in ascending alpha order for consistency
    wells.sort_index(inplace=True)
    names, data = [], []

    # define plot y axis max/min zscore by max value in data
    if isinstance(maxv, str) and maxv =='auto':
        maxval = 0
        for d in data:
            val1 = max(abs(d[0]))
            val2 = max(abs(d[1]))
            submax = max(val1, val2)
            if submax > maxval:
                maxval = submax
        maxv = round(maxval) + 2

    # create figure
    try:
        if len(wells.columns) == 0:
            print('no wells for skyline')
            return None
        if len(wells.columns) > 1:
            wells['consensus'] = pa.consensus(wells)
        l = len(wells.columns)
    except AttributeError:
        l = 1
        wells = pd.DataFrame(wells)

    fig = plt.figure(figsize=(10, (l*2)+2), facecolor="white")

    for i, col in enumerate(wells.columns):
        w = wells[col]
        # define subplot in grid, counting down from the max length
        ax = fig.add_subplot(l, 1, l-i, facecolor="white")
        ax.axis([0, 978, -maxv, maxv])
        ax = format_skyline_plot(ax)
        ax.set_ylabel(col)
        ax.set_yticks([-maxv, 0, maxv])
        # plot up + down values in green/red
        ax.bar(range(978), w.clip(lower=0), width=2, facecolor='green', edgecolor='green')
        ax.bar(range(978), w.clip(upper=0), width=2, facecolor='firebrick', edgecolor='firebrick')
        # plot optional highlight genes
        if highlights:
            # check if highlight is a string, if so tuck in list so looping works
            if isinstance(highlights, str):
                highlights = [highlights]
            for g in highlights:
                try:
                    vx = np.where(w.index.values==g)[0][0]
                except:
                    vx = np.where(w.index.values==g)[0]
                ax.bar(vx, maxv, width=2, facecolor='lightgray', edgecolor='lightgrey')
                ax.bar(vx, -maxv, width=2, facecolor='lightgray', edgecolor='lightgrey')
                ax.bar(vx, w.loc[g], width=2, facecolor='lime', edgecolor='lime')
                ax.bar(vx, w.loc[g], width=2, facecolor='red', edgecolor='red')
        # white line at zero
        ax.bar(range(978), np.zeros(978), width=1, facecolor='grey', edgecolor='grey')

    plt.suptitle(title)
    plt.savefig(os.path.join(outpath + '.png'), facecolor="white")
    plt.close()




def main():
    pass


if __name__ == '__main__':
    main()
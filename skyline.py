#!/usr/bin/python
'''
1/28/17 - separeted out from the growing plateanalysis script 
into justh the skyline plotting script, calling functions 
from plateanalysis as needed
'''
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import gct, os, pt, gt
import plateanalysis


def get_mergemode(f):
    """ look through map file info for a gct file to determine relationship amongst samples
    and replicates - within or across batches, with or without dose """
    g = gct.Gct(f)
    g.get_headers()
    h = g.fh
    numreps = []
    if 'dose' in h.columns:
        uniques = []
        for n in g.names:
            for d in h[h['name']==n]['dose'].unique():
                uniques.append((n, d))
        for b in g.batches:
            for u in uniques:
                n, d = u
                num = len(h[(h['batch']==b) & (h['name']==n) & (h['dose']==d)].index)
                numreps.append(num)
        if sum(numreps)/len(numreps) > 1.5:
            bmode = 'withinbatch'
        else:
            bmode = 'acrossbatch'
    else:
        for b in g.batches:
            for n in g.names:
                num = len(h[(h['batch']==b) & (h['name']==n)].index)
                numreps.append(num)
        if sum(numreps)/len(numreps) > 1.5:
            bmode = 'withinbatch'
        else:
            bmode = 'acrossbatch'
    if 'dose' in h.columns:
        if bmode == 'acrossbatch':
            mmode = 'plate_dose'
        elif bmode == 'withinbatch':
            mmode = 'batch_dose'
    else: 
        if bmode == 'acrossbatch':
            mmode = 'plate'
        elif bmode == 'withinbatch':
            mmode = 'batch'
    return mmode


def get_rep_wids(f, well, mmode):
    """ given a gct file, a particular well and the mergemode for that plate, return a
    list of well ids which are replicates of the input well """
    g = gct.Gct(f)
    g.get_headers()
    h = g.fh
    w = h.loc[well]
    if mmode == 'plate_dose':
        n = w['name']
        d = w['dose']
        wids = h[(h['name']==n) & (h['dose']==d)].index.values
    elif mmode == 'plate':
        n = w['name']
        wids = h[h['name']==n].index.values
    elif mmode == 'batch':
        n = w['name']
        b = w['batch']
        wids = h[(h['batch']==b) & (h['name']==n)].index.values
    elif mmode == 'batch_dose':
        n = w['name']
        b = w['batch']
        d = w['dose']
        wids = h[(h['batch']==b) & (h['name']==n) & (h['dose']==d)].index.values
    else:
        print('incorrect mmode passed')
        wids = []
    return wids


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


def plot_skylines(wells, highlights=None, outpath='dflt', t=1, maxv=10, title='auto', line1=None, line2=None):
    """ plot a series of individual skylines, ignoring any replicates/consensus. can be list of wells or
    dataframe, in which case each column (sample) is plotting individually """
    if isinstance(wells,list):
        for x in wells:
            try:
                plot_skyline(x, highlights=highlights, outpath=outpath, maxv=maxv)
            except KeyError:
                print(x, 'not found in index')
    else:
        for x in wells.columns:
            try:
                plot_skyline(wells[x],highlights=highlights, outpath=outpath, maxv=maxv)
            except KeyError:
                print(x, 'not found in index')


# add support for dictionary of highlights to have different classes/colors of highlight genes
# with a labeled legend
def plot_skyline(wells, highlights=None, outpath='dflt', t=1, maxv=10, title='auto', line1=None, line2=None):
    """ well can either be a single instance (pd Series) or a list of replicate instances, or dataframe,
     in which case the consensus skyline will be generated and plotted. maxv is yscale, 10 dflt or 'auto'
     which sets based upon max value in the dataset. line 1 and line 2 may contain metadata
     for labelling the plot, and highlights accepts a list of genes to be accented """
    # check if wells is long its just a single series, and put it in a list container
    # otherwise if wells is multiple replicate instances (length 2-5) leave as is
    try:
        if len(wells.columns) > 1:
            wells = [wells.loc[:,i] for i in wells.columns]
    except AttributeError:
        if len(wells) > 10:
            wells = [wells]
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
    for w in wells:
        data.append(sep_updown(w, t=t))
        names.append(w.name)
    # pass to external function to get consensus of wells if multiple instances
    if len(wells) > 1:
        c = well_consensus(wells, 'consensus')
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
    if outpath is 'dflt':
        outpath = gt.dflt_outpath(fldr_name='output figs')
        if title is 'auto':
            fn = wells.iloc[:,0].name.replace(':','-')
        else:
            fn = title
        outpath = os.path.join(outpath, fn)
    # sort genes in ascending alpha order for consistency
    wells.sort_index(inplace=True)
    names, data = [], []

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
    try:
        wells['consensus'] = gt.consensus(wells)
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

    # load in metadata from optional line arguments and incorporate into title on multiple lines
    if title is 'auto':
        title = pt.get_shn(outpath)
    if info is not None:
        title = title + '\n' + str(info)

    plt.suptitle(title)
    plt.savefig(outpath, facecolor="white")
    plt.close()


# should allow a list not a file to be passed in -- Feb 16 suggestion
def selective_skyline(path, opath, welllist, highlights):
    """ path as gct directory, outpath for figure files, well list is text file in comma delim format
    which permits same or different plates PGA314:A14, PGA318:G09
    could use a script to generate that list of all replicates from a map file """
    flist = pt.get_flist(path, '.gct')
    shnames = {}
    for f in flist:
        shnames[f.split('/')[-1].split('_',1)[0]] = f
    prev_plates = ''
    with open(welllist, 'rU') as infile:
        for line in infile:
            plates = list(set([p.split(':')[0] for p in line.strip().split(',')]))
            wells = line.strip().split(',')
            if plates != prev_plates:
                if len(plates) > 1:
                    files = [shnames[n] for n in plates]
                    gd, dfd = gct.load_panel(files, indv=True)
                    plist = []
                    for d in dfd.values():
                        plist.append(gt.addr_id_df(d))
                    df = pd.concat(plist, axis=1)
                elif len(plates) == 1:
                    g = gct.Gct(shnames[plates[0]])
                    df = g.build_dframe()
                    df = gt.addr_id_df(df)
            prev_plates = plates
            wids = []
            outpath = os.path.join(opath, wells[0].replace(':', '-'))
            for w in wells:
                try:
                    wids.append(df[w])
                except KeyError:
                    print(w, ' not in columns, len:', len(df.columns))
            plot_skyline(wids, outpath, highlights)


def selective_skyline_meta(path, opath, welllist, highlights):
    """ input file should come from 'annotate wells' with appropriate metadata fields of form:
    name    dose    rank    zs    cc prom    PGA103:A14, PGA104:A14 (one or many wells ok) """
    # get files with associated short plate name dict, remove rank files
    flist, shnames = pt.get_flist(path, '.gct', shn=True)
    flist = [x for x in flist if 'rank' not in x]
    prev_plates = ''
    with open(welllist, 'rU') as infile:
        for line in infile:
            line = line.strip().split('\t')
            wells = line[-1].split(',')
            wells = [w.strip() for w in wells]
            metad = line[:-1]
            plates = list(set([p.split(':')[0] for p in wells]))
            hd = ['name:', 'dose:', 'rank: ', 'zs: ', 'cc: ', 'prom: ']
            # associate metadata fields present in above order to headers
            line = list(zip(hd, metad))
            # clean up string formatting from python list output
            # split line into the two different header lines for the skyline plot
            line1, line2 = line[:2], line[2:]
            line2 = ['{:.2f}'.format(float(x)) if '.' in x else x for x in line2]
            line1 = str(line1).replace("'", "")
            line1 = line1.replace(',', '')
            line2 = str(line2).replace("'", "")
            line2 = line2.replace(',', '')
            # check if same set of plates in use, otherwise load file or panel
            if plates != prev_plates:
                if len(plates) > 1:
                    files = [shnames[n] for n in plates]
                    gd, dfd = gct.load_panel(files, indv=True)
                    plist = []
                    # assemble consensus dataframe
                    for d in dfd.values():
                        plist.append(gt.addr_id_df(d))
                    # concat all wells
                    df = pd.concat(plist, axis=1)
                elif len(plates) == 1:
                    g = gct.Gct(shnames[plates[0]])
                    df = g.build_dframe()
                    df = gt.addr_id_df(df)
            # re-define previous plates
            prev_plates = plates
            wids = []
            # define file destination
            # changed to include name and dose
            name = line[0][1]
            if name.startswith('na') or name == '-666':
                name = wells[0]
            dose = line[1][1]
            if dose != 'na':
                dose = '{:.2f}'.format(float(dose) * 1E6)
            desc = plates[0] + '-' + name + '-' + str(dose)
            print(desc)
            # if the name is just the well, look for colon and use that
            # otherwise use the name + dose
            if ':' in desc:
                outpath = os.path.join(opath, wells[0].replace(':', '-'))
            else:
                outpath = os.path.join(opath, desc.replace(':', '-'))
            outpath = outpath.replace('.', ',') + '.png'
            print(outpath)
            for w in wells:
                try:
                    wids.append(df[w])
                except KeyError:
                    print(w, ' not in columns, len:', len(df.columns))
            args = [wids, outpath, line1, line2, highlights]
            plot_skyline(*args)


def batch_run_meta_skyline():
    """ paste values here and run this to do so in batch mode easily """
    path = '/Volumes/WRBHDD/wrb/Desktop/all processed/'
    welllist = '/Volumes/WRBHDD/wrb/Desktop/MKD303wells.txt'
    platelist = '/Volumes/WRBHDD/wrb/Desktop/overall mergelist (8 June 2017).txt'
    gene = '200678_x_at'
    opath = '/Volumes/WRBHDD/wrb/Desktop/skylines/'
    outpath = welllist.replace('.txt', '-reps.txt')

    # plateanalysis.clean_zsfiles(path)

    plateanalysis.get_reps(path, welllist, platelist, outpath)

    welllist = outpath
    plateanalysis.annotate_wells(path, welllist, gene)

    welllist = welllist.replace('.txt', '-annot.txt')
    skyline.selective_skyline_meta(path, opath, welllist, gene)


def main():
    pass


if __name__ == '__main__':
    main()
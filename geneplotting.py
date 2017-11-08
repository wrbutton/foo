#!/usr/bin/python
"""
12/10/16 - plotting script version 3.0

multipurpose plotting script to plot genes from gct files in
a variety of configurations specified from the command line,
and described using the argparse help.

gct files of interest and a genelist.txt file containing genes

are placed in the input folder (Desktop/plotfiles by default)
and images created and written to output folder.

options are provided to plot by dotplot or plate view, with panels
of etiher plates or genes, with configurable options

terminology:
vctr = data vector, one row of gene values
dfl = Gct class object of data file, including header info etc
dfr = pandas data frame created from gct, based on gct
spl = subplot, each component figure

refer to argparse help for descriptions of options and commands
"""
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys, os, string, argparse
import gct, gt


def get_files(filepath, ext):
    # build list of files in targeted directory with specified extension
    file_list, flist = [], []
    for input_file in os.listdir(filepath):
        fullfilename = os.path.abspath(os.path.join(filepath, input_file))
        if input_file.endswith(ext):
            file_list.append(fullfilename)
    return file_list


def get_awells():
    # returns list of 384 three character well IDs
    awells = []
    rows = string.ascii_uppercase[0:16]
    cols = range(1,25)
    for l in rows:
        for n in cols:
            # join lettrs and nums with 2 characater num format
            awells.append(l + str('{:02d}'.format(n)))
    return awells


def txt2list(genefile):
    # assemble a list from text file with one item per line
    mylist = []
    with open(genefile, 'rU') as in_file:
        for line in in_file:
            mylist.append(line.strip())
    return mylist


def prep_file(file):
    # use gct class to get header info and pd.dataframe
    dfl = gct.Gct(file)
    dfr, h = dfl.build_dframe()
    return dfl, dfr


def prep_vctr(dfr, feature, order):
    # grab the data row for given feature, and align to full 384-well length
    vctr = dfr.loc[feature]
    # vctr.index = vctr.index.map(lambda x: x.split(':')[1])
    new_index = get_awells()
    if any(vctr.index.str.contains(':')):
        pname = vctr.index.values[0].split(':')[0]
        new_index = [pname + ':' + x for x in new_index]
    awells = pd.Series(np.nan * 384, index=new_index)
    # merge the new wells onto an array of 384 in order to keep spacing
    fvctr = awells.combine_first(vctr)
    if order == 'col':
        print('plotting in column order')
        # create a new sorted vector by sorting by number first, then row
        svctr = fvctr.loc[sorted(fvctr.index.values, key=lambda x: (x[1:],x[0]))]
    else:
        # sort vector samples by well row value by default
        svctr = fvctr.loc[sorted(fvctr.index)]
    return svctr


def save_fig(fig, filename, figsize):
    #fig.set_tight_layout(True)
    fig.set_size_inches(figsize[0], figsize[1])
    fig.savefig(filename)
    plt.close()


def format_axes(myaxes):
    # formats the suplot for dotplot, removing figure boundaries
    # myaxes.grid(True)
    myaxes.set_facecolor('white')
    myaxes.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    myaxes.tick_params(axis='y', which='both', left='on', right='off', labelleft='on')
    for item in myaxes.get_yticklabels():
        item.set_fontsize(8)
    # leave only the left y axis vertical line, remove others
    myaxes.spines['left'].set_visible(True)
    myaxes.spines['right'].set_visible(False)
    myaxes.spines['top'].set_visible(False)
    myaxes.spines['bottom'].set_visible(False)


def re_order(mylist, orderlist):
    # re order gene or file list to desired order
    if len(set(orderlist.split(','))) != len(mylist):
        sys.exit('oops, order list error! use comma separated integers')
    norderlist = [int(n) - 1 for n in orderlist.split(',')]
    a = np.array(mylist)
    newlist = list(a[norderlist])
    return newlist


def plot_multigene(pdim, ptype, file_list, genelist, outpath, figsize, order='row', stype=True, test=False):
    print('plotting multiple genes per plate with', ptype, 'plot')
    if order == 'col':
        print('plotting in column order')
    for file in file_list:
        # open one file at a time, pull genes and plot subplots
        dfl, dfr = prep_file(file)
        fig = plt.figure(frameon=False)
        fig.suptitle(dfl.shortname, fontsize=16)
        for idx, feature in enumerate(genelist):
            # increment through subplots per gene
            spl = plt.subplot(pdim.r, pdim.c, idx+1)
            print('subplot = ', pdim.r, pdim.c, idx+1)
            vctr = prep_vctr(dfr, feature, order)
            if ptype == 'dot':
                make_dotplot(dfl, spl, feature, vctr, stype)
                suffix = '_dot.png'
            elif ptype == 'plate':
                make_plateplot(spl, feature, vctr)
                suffix = '_plate.png'
        filename = os.path.join(outpath, feature + suffix)
        plt.tight_layout()
        save_fig(fig, filename, figsize)
        foo
        if test == True:
            sys.exit('test mode, quitting after first figure')


def plot_multiplate(pdim, ptype, file_list, genelist, outpath, figsize,
                                        order='row', stype=True, test=False):
    print('plotting multiple plates per gene with', ptype, 'plot')
    if order == 'col':
        print('plotting in column order')
    data_list = []
    # open and load all files, keep dframes in memory
    for file in file_list:
        dfl, dfr = prep_file(file)
        data_list.append((dfl, dfr))
    # then work through gene list, pulling from all plates
    for feature in genelist:
        fig = plt.figure(frameon=False)
        fig.suptitle(feature, fontsize=16)
        for idx, (dfl, dfr) in enumerate(data_list):
            spl = plt.subplot(pdim.r, pdim.c, idx+1)
            vctr = prep_vctr(dfr, feature, order)
            if ptype == 'dot':
                make_dotplot(dfl, spl, dfl.shortname, vctr, stype)
                suffix = '_dot.png'
            if ptype == 'plate':
                make_plateplot(spl, dfl.shortname, vctr)
                suffix = '_plate.png'
        filename = os.path.join(outpath, feature + suffix)
        save_fig(fig, filename, figsize)
        if test == True:
            sys.exit('test mode, quitting after first figure')


def make_plateplot(spl, splabel, vctr):
    # set additional title and axis
    plt.ylabel(splabel, fontsize=16)
    row_labels = list(string.ascii_uppercase[0:16])
    row_range = list(np.arange(16))
    plt.yticks(row_range, row_labels, fontsize=8)
    col_labels = list(np.arange(1,25))
    col_range = list(np.arange(0,24))
    plt.xticks(col_range, col_labels, fontsize=9)
    plt.tick_params(labelright=True, labeltop=True)
    # this sets the tick length to zero, but leaves labels
    plt.tick_params(axis=u'both', which=u'both',length=0)
    # reshape array and plot
    d = vctr.values.reshape(16,24)
    im = plt.imshow(d, interpolation='nearest', cmap='inferno')
    # use matplotlib axes1 to keep colorbars in line with figs
    divider = make_axes_locatable(spl)
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


def make_dotplot(dfl, spl, splabel, vctr, stype):
    format_axes(spl)
    # set additional title and axis
    plt.ylabel(splabel, fontsize=12)
    if vctr.min() > 0:
        plt.axis([-2, 390, 2, 16])
    elif vctr.min() < 0:
        if vctr.max() > 10 or vctr.min() < -10:
            plt.axis([-2, 390, -20, 20])
        else:
            plt.axis([-2, 390, -10, 10])
    awells = list(vctr.index)
    # plot primary data, and then poscons
    if stype == False:
        plt.plot(vctr.values, color='blue', marker='o', ls='', markersize=5, mew=0)
    elif stype == True:
        vehy = vctr[dfl.vehicles]
        vehx = [awells.index(w) for w in dfl.vehicles]
        posy = vctr[dfl.poscons]
        posx = [awells.index(w) for w in dfl.poscons]
        # foobar
        plt.plot(vctr.values, color='silver', marker='o', ls='', markersize=5, mew=0)
        plt.plot(vehx, vehy, color='b', marker='o', ls='', markersize=5, mew=0)
        plt.plot(posx, posy, color='r', marker='o', ls='', markersize=5, mew=0)
    elif stype == 'other':
        searchwells = dfl.names['unstimulated']
        testy = vctr[dfl.test]
        testx = [awells.index(w) for w in dfl.test]
        plt.plot(testx, testy, color='silver', marker='o', ls='', markersize=5, mew=0)
        vehy = vctr[dfl.vehicles]
        vehx = [awells.index(w) for w in dfl.vehicles]
        plt.plot(vehx, vehy, color='b', marker='o', ls='', markersize=5, mew=0)
        othery = vctr[searchwells]
        otherx = [awells.index(w) for w in searchwells]
        plt.plot(otherx, othery, color='g', marker='o', ls='', markersize=5, mew=0)


def main():
    print('') # spacer line
    parser = argparse.ArgumentParser()
    # required arguments / choices first
    parser.add_argument('dims', nargs=2, type=int,
                 help='subplot grid dimensions, enter: rows <space> columns')
    group1 = parser.add_argument_group('plot layout')
    group1.add_argument('-mp', '--multiplate', action='store_true',
                 help='one plot per gene, multiple plate subplots')
    group1.add_argument('-mg', '--multigene', action='store_true',
                 help='one plot per plate, multiple gene subplots')
    group2 = parser.add_argument_group('plot format')
    group2.add_argument('-dp', '--dotplot', action='store_true',
                 help='horizontal dotoplot of gene expression')
    group2.add_argument('-pp', '--plateplot', action='store_true',
                 help='plate based plot of gene expression')
    # optional flags and definitions
    parser.add_argument('-c', '--columnorder', action='store_true',
                 help='plot in column order, row dflt')
    parser.add_argument('-t', '--test', action='store_true',
                 help='test mode, only generates one image to check')
    parser.add_argument('-ds', '--dontsep', action='store_true',
                 help="don't separate vehicles and poscons, plot one series")
    parser.add_argument('-os', '--othersep', action='store_true',
                 help="dont plot poscons, highlight vehicles and provided test name")
    parser.add_argument('-sl', '--speclayout', action='store_true',
                 help='specify layout of suplot images, start at 1, no blanks')
    parser.add_argument('-fl', '--filelist', action='store_true',
                 help='''use data files in /plotfiles/filelist.txt
                 with the text file listing absoltue paths of data files''')
    parser.add_argument('--fsize', nargs=2, type=int, default=(15,8),
                 help='size in inches for final fig, eg: 15 8')
    group2 = parser.add_argument_group('locations')
    group2.add_argument('--input', help='''input data file dir including 
                 genefile.txt, default: Desktop/plotfiles/''',
                 default=gt.dflt_outpath(fldr_name='plotfiles'))
    group2.add_argument('--output', help='''output directory for images, 
                 default: Desktop/output_images/''',
                 default=gt.dflt_outpath(fldr_name='output_images'))
    args = parser.parse_args()

    # create output folder if it doesn't exist
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    # assemble optional keyword arg dictionary to pass
    kwargs = {}
    if args.columnorder:
        kwargs['order'] = 'col'
    if args.dontsep:
        kwargs['stype'] = False
    if args.test:
        kwargs['test'] = True
    if args.othersep:
        kwargs['stype'] = 'other'

    # assemble basic argumnets of locations and lists
    filepath, outpath = args.input, args.output
    genefile = os.path.join(filepath, 'genelist.txt')
    if not os.path.exists(genefile):
        sys.exit('error: expected genelist.txt not present in input folder ' + filepath)
    genelist = txt2list(genefile)

    # set the list of data files, either specified list or in folder
    if args.filelist:
        file_list = txt2list(os.path.join(filepath, 'filelist.txt'))
    else:
        file_list = get_files(filepath, '.gct')

    # print out genelist and filelists in use, except if long
    fnames = [f.split('/')[-1] for f in file_list]
    if len(fnames) <= 16:
        print('file list: ', fnames)
    else:
        print('file list longer than 16 files')
    if len(genelist) <= 16:
        print('gene list: ', genelist)
    else:
        print('gene list longer than 16 genes')
    print('subplot dimensions: ', args.dims)

    # define named tuples for more convenient later access
    Dim = namedtuple('Dim', 'r, c')
    pdim = Dim(args.dims[0], args.dims[1])
    Size = namedtuple('Size', 'w, h')
    figsize = Size(args.fsize[0], args.fsize[1])

    # determine plotting format
    if args.dotplot:
        ptype = 'dot'
    elif args.plateplot:
        ptype = 'plate'
    else:
        print('select a plot format: dotplot(dp) or plateplot(pp)')

    # set up reorderd list for specified layouts
    if args.speclayout:
        print('\n original order:')
        if args.multigene:
            for g in genelist:
                print(g)
            print('input comma sep list of genes in desired order')
            print('eg. A, B, C w/ reorder list 3, 1, 2 = C, A, B')
            neworder = input('re-order list: ')
            genelist = re_order(genelist, neworder)
        elif args.multiplate:
            for p in file_list:
                print(p)
            print('input comma sep list of plates in desired order')
            print('eg. A, B, C w/ reorder list 3, 1, 2 = C, A, B')
            neworder = input('re-order list: ')
            file_list    = re_order(file_list, neworder)

    # pack argument list to pass to function call
    argmnts = (pdim, ptype, file_list, genelist, outpath, figsize)

    # determine plotting layout and call functions to make plots
    if args.multigene:
        plot_multigene(*argmnts, **kwargs)
    elif args.multiplate:
        plot_multiplate(*argmnts, **kwargs)
    else:
        print('select plot layout: multigene(mg) or multiplate(mp)')

if __name__ == '__main__':
    main()

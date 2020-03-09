#!/usr/bin/python

""" higher level functions tying together commands from other modules to consolidate the 
streamlined  processing of basic dataset analysis and plotting """

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys, os, glob, gct, gcsv
import qc, imgs
import warnings
warnings.filterwarnings(action='once')
import gt, pa, dim_reduct
import plottools, pickle, sys
import collections as cll





def run_plate_analysis(mode='ind', cats='nd', path='dflt'):
    """ runs standard analysis on either each plate individually 'ind' or all togegther 'comb'
    most useful for plates with doses. the default loc  

    default path will be newQC on the desktop """

    path = gt.check_dfltarg(path, os.path.join(gt.check_desktop(), 'newQC'))

    fl = gt.globit(path, '*ZSVCQNORM*')

    print(fl)

    if mode == 'comb':
        dl, hl = [], []
        for i,f in enumerate(fl):
            d, h = gct.extractgct(f)
            if i == 0:
                try:
                    pname = d.name + '+'
                except:
                    pname = h.addr[0].split(':')[0] + '+'
            if len(h.batch.unique()) > 1:
                # fix sample labels for plate/batch
                h.plate = h.plate + h.batch
            # define labels (should I add plate?)
            h = gt.gen_label(h, cats)
            dl.append(d)
            hl.append(h)
        try:
            d = pd.concat(dl, axis=1)
            d.name = pname
        except ValueError:
            sys.exit('no gct file plates to analyze')
        h = pd.concat(hl, axis=0)

        analyze_plate(d, h, cats)

    elif mode == 'ind':
        for f in fl:
            d, h = gct.extractgct(f)
            # define labels (should I add plate?)
            h = gt.gen_label(h, cats)

            analyze_plate(d, h, cats)



def analyze_plate(d, h, cats):
    """ take a dataset (assuming with doses) and generate standard output figs """
    # should worry about organizing output too

    # create consensus, create SC plot
    dc, hc = pa.assemble_consensus(d, h, cats, save=True, sc=True)
    hc = gt.gen_label(hc, 'nd')
    dc.name = d.name

    # tSNE simple first pass, two parameters
    dim_reduct.tsne2(dc, hc, px=10, lr=[10, 150], inter=True)

    # create general dendrogram (only, no heatmap)
    dim_reduct.make_dendrogram(dc, labels=hc.label, outpath=True)

    # plot correlation matrix of the sorted combined zs by name and dose
    # can then follow up to plot the sweep or the clustered
    plottools.plot_correlation_matrix(dc, hc, title='dflt', sort=True, outpath=True, sparselabel=True, grid=True, labels=hc.name)

    if 'd' in cats:
        # plot landmark concs
        newcats = cats.replace('d','')
        outpath = gt.dflt_outpath(fldr_name='landmark concs')
        #plottools.plot_landmark_concs(dc, hc, cats=newcats, genes='all', labels='dose', outpath=outpath)

        # call combo function to find genes that move and plot those dose-response plots (30 per plate)
        plottools.plot_ex_genes(d, h, n=10, mode='med')



#!/usr/bin/python

""" functions to assemble data across plates and sets into larger denominations """

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gct, gt, os, gt



def gather_ref(path):
    """ assemble reference wells across plates and save as file """
    fl = glob.glob(path + '/**/*_ref_*.gct', recursive=True)
    datlist = []
    for f in fl:
        d, h = gct.extractgct(f)
        newcols = [x for x in d.columns if 'B01' not in x]
        nd = d[newcols].mean(axis=1)
        nd.name = d.columns[0].split(':')[0]
        datlist.append(nd)
    data = pd.concat(datlist, axis=1)
    outpath = os.path.join(path, 'ref_consolidation.txt')
    data.to_csv(outpath, sep='\t')


def consolidate_vehicles(path):
    """ assemble reference wells across plates and return dataset """
    sstring = '*QNORM_sorted*.gct'
    data = []
    fl = gt.globit(path, sstring)
    for f in fl:
        d, h = gct.extractgct(f)
        ds = gt.dosub(d, h, {'type':'vehicle'}).sort_index()
        ds = ds.mean(axis=1)
        ds.name = d.name
        data.append(ds)
    fdata = pd.concat(data, axis=1)
    return fdata


def assemble_poscons_from_files(path):
    """ assemble reference wells across plates and return dataset """
    data = pd.DataFrame()

    flist = gt.get_flist(path, '.gct')

    for f in flist:
        d, h = gct.extractgct(f)
        sh = h[h['type']=='poscon']
        wids = sh.index.values
        labels = sh['name'].values
        sub = d[wids]
        sub.columns = [x + ':' + y for x, y in zip(sub.columns, labels)]
        data = pd.concat([data, sub], axis=1)

    # data.to_excel(os.path.join(path, 'data.xlsx'))
    return data


def assemble_data_from_files(path, arg_dict):
    """ takes in dictionary of {param: value} to scrape the given gct file and return data meeting ALL"""
    d, h = gct.extractgct(path)
    sh = h
    for c, v in arg_dict.items():
        if isinstance(v, str):
            sh = sh[sh[c]==v]
        elif isinstance(v, list):
            sh = sh[sh[c].isin(v)]
        else:
            print('error with filter dictionary value')
    if len(sh.index.values) == 0:
        print('no wells in selection')
    print(len(sh.index.values))
    return d[sh.index.values], sh


def cheap_assemble_data_from_files(path, arg_dict):
    """ takes in dictionary of {param: value} to scrape the given gct file and return data meeting ALL"""
    g = gct.Gct(path)
    g.get_headers()
    sh = gt.hsub(g.fh, arg_dict)
    ds, hs = g.build_sample_subset(sh.index.values)
    print(str(len(sh.index.values)) + ' wells taken from ' + gt.get_shn(path))
    return ds, sh


def grab_data_from_folder(path, arg_dict, save=False, outpath='dflt'):
    """ takes in dictionary of {param: value} to scrape the given folder for all gct files and return data meeting ALL"""
    flist = gt.get_flist(path, '.gct')
    dlist, hlist = [], []
    for f in flist:
        sd, sh = cheap_assemble_data_from_files(f, arg_dict)
        dlist.append(sd)
        hlist.append(sh)
    d = pd.concat(dlist, axis=1)
    h = pd.concat(hlist, axis=0)
    d.index.name = ''
    if save is True:
        fold = os.path.basename(os.path.dirname(path))
        name = ''
        for k, v in arg_dict.items():
            name += ',{}-{}'.format(k,v)
        name += '.gct'
        name = fold + name
        print(name)
        outpath = gt.dflt_outpath(outpath, 'foo', fn=name)
        print(outpath)
        gct.save_headergct(d, outpath, h)
    else:
        return d, h


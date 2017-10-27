#!/usr/bin/python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gct, pt, os



def assemble_poscons_from_files(path):
    data = pd.DataFrame()

    flist = pt.get_flist(path, '.gct')

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
    arglist = []
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



def grab_data_from_folder(path, arg_dict):
    """ takes in dictionary of {param: value} to scrape the given folder for all gct files and return data meeting ALL"""
    flist = pt.get_flist(path, '.gct')
    dlist, hlist = [], []
    for f in flist:
        sd, sh = assemble_data_from_files(f, arg_dict)
        dlist.append(sd)
        hlist.append(sh)
    d = pd.concat(dlist, axis=1)
    h = pd.concat(hlist, axis=0)
    d.index.name = ''
    # if rename == True:
    #     h['hdr'] = h.index + ':' + h.name
    #     d.columns = h.hdr
    # elif rename == 'dose':
    #     h['hdr'] = h.index + ':' + h.dose + ':' + h.name
    #     d.columns = h.hdr
    return d, h


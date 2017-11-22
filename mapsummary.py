#!/usr/bin/python -tt

import gt, os, plottools
import collections as cll
import pandas as pd
import numpy as np
import gct


def ctup(df, arg_dict, col=':', u=True, lst=False):
    if lst is True and u is True:
        r = gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna().unique()
    elif lst is False and u is True:
        r = len(gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna().unique())
    elif lst is True and u is False:
        r = gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna()
    elif lst is False and u is False:
        r = len(gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna())
    return r

    

def plate_map_vis(myseries):
    """ just translate directly into a dict or overwrite to use?
    returns just the array to plot """
    xltr = dict(zip(myseries.unique(), range(0,len(myseries.unique()))))
    data = myseries.dropna().apply(lambda x: xltr[x])
    others = []
    data = data.reindex(data.index.values + others)
    plottools.plot_plateplot(data, outpath='dflt', label=data.name, cmap='tab20', clrbar=xltr)


def check_maps(path):

    pert_dict, map_list = {}, {}
    wellpert_dict = {}



    flist = gt.get_flist(path, ext='.xlsx')
    flist = [x for x in flist if len(os.path.split(x)[-1])==11]
    awells = gt.get_awells()

    headers = {'wells #': lambda x: len(x.index),
               'test #': lambda x: ctup(x, {'type':'test'}, 'well'),
               'vehicle #': lambda x: ctup(x, {'type':'vehicle'}, 'well'),
               'poscon #': lambda x: ctup(x, {'type':'poscon'}, 'well'),
               'dose #': lambda x: ctup(x,{'type':'test'},'dose [M]',u=False)/ctup(x, {'type':'test'}, 'name', u=False),
               'poscons': lambda x: ctup(x, {'type':'poscon'}, 'name', lst=True)}

    summary = pd.DataFrame(columns=headers)

    for f in flist:
        pname = gt.get_shn(f).split('.')[0]
        print(pname)
        m = pd.read_excel(f)
        m.sort_index(inplace=True)
        batches = m['batch'].dropna().unique()
        # check wells for full plate
        if len(set(awells) - set(m['well'].values)) != 0:
            print('{} wells error, {} entries'.format(pname, len(m.index)))
        # check vehicle named 'DMSO'
        veh_names = m.loc[m['type']=='vehicle','name'].unique()
        if veh_names != 'DMSO':
            print('{} vehicle error, names: {}'.format(pname, veh_names))
   


        # summarize the header info per batch, and assemble pert-lists
        # for the overlap comparisons
        for b in batches:
            entry = pname + '-' + b
            ms = gt.hsub(m, {'batch': b})
            # gather pert names for overlap comparison
            pert_dict[entry] = ctup(m, {'batch':b, 'type':'test'}, 'name', lst=True)
            # get the well-pert identities for same plate comparison
            ms['addr'] = ms['well'] + '-' + ms['name'].apply(lambda x: str(x))
            wellpert_dict[entry] = ms['addr'].values
            for k in headers.keys():
                summary.loc[entry, k] = headers[k](ms)




    
    summary.to_excel(os.path.join(path, 'batch_composition.xlsx'))
    same_plates = gt.overlap_matrix(wellpert_dict.values(), wellpert_dict.keys())
    pert_overlap = gt.overlap_matrix(pert_dict.values(), pert_dict.keys())
    pert_overlap.to_excel(os.path.join(path, 'pert_overlaps.xlsx'))
    same_plates.to_excel(os.path.join(path, 'well-name_overlaps.xlsx'))
    
    """
    create a map viz template for each plate (after wells checked) 
    maybe helpful:
    dc = m[m.columns[m.columns.str.contains('dose')==True]].columns
    """


def main():
    path = '/Users/wrb/Desktop/AVA308_ZSVCQNORM_n358x978.gct'
    h = gct.extractheader(path)
    plate_map_vis(h['type'])


if __name__ == '__main__':
  main()

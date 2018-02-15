#!/usr/bin/python -tt

import gt, os, plottools
import collections as cll
import pandas as pd
import numpy as np
import gct, math


def compare_files(path1, path2, show=False):
    difflines = []
    line_count = 0
    shn1, shn2 = gt.get_shn(path1), gt.get_shn(path2)
    with open(path1, 'r') as f1, open(path2, 'r') as f2:
        for l1 in f1:
            l2 = f2.readline()
            if l1 == l2:
                continue
            else:
                line_count += 1
                difflines.append(l1)
    print(shn1, shn2)
    if show is True:
        print(' '.join(difflines))
    else:
        print(line_count)


def compare_samename_files(path1, path2, show=False):
    fl1, fl2 = gt.get_flist(path1), gt.get_flist(path2)
    for file in fl1:
        shn = gt.get_shn(file)
        other = os.path.join(path2, shn)
        if os.path.isfile(other):
            compare_files(file, other, show=show)


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


def plate_map_vis(myseries, path='dflt'):
    """ just translate directly into a dict or overwrite to use?
    returns just the array to plot """
    xltr = dict(zip(myseries.unique(), range(0,len(myseries.unique()))))
    #data = myseries.dropna().apply(lambda x: xltr[x])
    data = myseries.apply(lambda x: xltr[x])
    if path is 'dflt':
        outpath = gt.dflt_outpath(fn=myseries.name)
    else:
        outpath = os.path.join(path, myseries.name)
    plottools.plot_plateplot(data, outpath=outpath, label=data.name, cmap='tab10', clrbar=xltr)


def check_maps(path):

    plot_fields = ['type', 'batch']
    # add checks to add batch and dose if present

    pert_dict, map_list = {}, {}
    wellpert_dict = {}

    flist = gt.get_flist(path, ext='.xlsx')
    flist = [x for x in flist if len(os.path.split(x)[-1]) == 11]
    if len(flist) == 0:
        flist = gt.get_flist(path, ext='.txt')
        flist = [x for x in flist if len(os.path.split(x)[-1])==10]
    print('flist = ', flist)

    dose_field = 'dose (M)'

    awells = gt.get_awells()

    headers = {'wells #': lambda x: len(x.index),
               'test #': lambda x: ctup(x, {'type':'test'}, 'well'),
               'vehicle #': lambda x: ctup(x, {'type':'vehicle'}, 'well'),
               'poscon #': lambda x: ctup(x, {'type':'poscon'}, 'well'),
               'dose #': lambda x: ctup(x,{'type':'test'},dose_field,u=False)/ctup(x, {'type':'test'}, 'name', u=False),
               'poscons': lambda x: ctup(x, {'type':'poscon'}, 'name', lst=True)}

    summary = pd.DataFrame(columns=headers)

    for f in flist:
        pname = gt.get_shn(f).split('.')[0]
        print(pname)
        if f.endswith('.xlsx'):
            m = pd.read_excel(f)
        elif f.endswith('.txt'):
            m = pd.read_table(f, index_col=False)
        m.sort_index(inplace=True)
        batches = m['batch'].dropna().unique()
        # check wells for full plate
        well_result = set(awells) - set(m['well'].values)
        if len(well_result) != 0:
            print('{} wells error, {} entries - {}'.format(pname, len(m.index), well_result))
        # check vehicle named 'DMSO'
        veh_names = m.loc[m['type']=='vehicle','name'].unique()
        #if len(veh_names) > 1 or not np.isnan(veh_names[0]):
        #    print('vehicles named: ', veh_names)
        #if veh_names != 'DMSO':
        #    print('{} vehicle error, names: {}'.format(pname, veh_names))

        # summarize the header info per batch, and assemble pert-lists
        # for the overlap comparisons
        for b in batches:
            entry = pname + '-' + b
            ms = gt.hsub(m, {'batch': b})
            # gather pert names for overlap comparison
            pert_dict[entry] = ctup(m, {'batch':b, 'type':'test'}, 'name', lst=True)
            # get the well-pert identities for same plate comparison
            ms.loc[:,'addr'] = ms['well'] + '-' + ms['name'].apply(lambda x: str(x))
            wellpert_dict[entry] = ms['addr'].values
            for k in headers.keys():
                try:
                    summary.loc[entry, k] = headers[k](ms)
                except KeyError:
                    if k is 'dose #':
                        dose_field = 'dose [M]'
                        try:
                            summary.loc[entry, k] = headers[k](ms)
                        except KeyError:
                            dose_field = 'dose'
                            try:
                                summary.loc[entry, k] = headers[k](ms)
                            except KeyError:
                                summary.loc[entry, k] = 'na'
                    else:
                        summary.loc[entry, k] = 'na'

        for pf in plot_fields:
            plot_series = m[pf]
            if len(plot_series.dropna().unique()) > 1:
                plot_series.name = pname + ' ' + pf
                plate_map_vis(plot_series, path=path)

    
    summary.to_excel(os.path.join(path, 'batch_composition.xlsx'))
    same_plates = gt.overlap_matrix(wellpert_dict.values(), wellpert_dict.keys())
    pert_overlap = gt.overlap_matrix(pert_dict.values(), pert_dict.keys())
    pert_overlap.to_excel(os.path.join(path, 'pert_overlaps.xlsx'))
    same_plates.to_excel(os.path.join(path, 'well-name_overlaps.xlsx'))



def main():
    path = '/Users/wrb/Desktop/AVA308_ZSVCQNORM_n358x978.gct'
    h = gct.extractheader(path)
    plate_map_vis(h['type'])


if __name__ == '__main__':
  main()

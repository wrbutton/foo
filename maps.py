#!/usr/bin/python -tt

import gt, os, plottools, string
import collections as cll
import pandas as pd
import numpy as np
import gct, math


def calc384well(well, batch):
    """ given a 96 well position in 3char form and the batch designation
    return the destination 384 well location
    A   B
    C   D   """

    # bump the index of alphabet up by one for the row references
    alpha = '_' + string.ascii_uppercase
    row, col = well[:1], int(well[1:])
    row = alpha.index(row)

    # sets up dictionaries for row and column adjusters
    bdict = dict(zip(alpha[1:5], [(-1, -1),
                                 (-1, 0),
                                 (0, -1),
                                 (0, 0)]))

    radj, cadj = bdict[batch]
    newrow = 2 * row + radj
    newcol = 2 * col + cadj
    if newcol < 0:
        newcol = 0
    if newrow < 0:
        newrow = 0
    newwell = alpha[newrow] + f'{newcol:02.0f}'
    return newwell


def merge96to384(flist):
    """pass in four txt/xlsx files named PCA102_A, PCA102_B... outputs full one in that dir
    layout: A   B
            C   D """
    if isinstance(flist, str):
        if '.txt' in flist:
            ext = '.txt'
        elif '.xlsx' in flist:
            ext = '.xlsx'
        else:
            print('map extension error')
        flist = gt.splitpaths(flist, ext=ext)

    fullplate = pd.DataFrame()

    namelist = []

    for b in string.ascii_uppercase[:4]:
        try:
            file = [x for x in flist if f'_{b}' in x][0]

            if ext == '.xlsx':
                m = pd.read_excel(file)
            elif ext == '.txt':
                m = pd.read_excel(file, sep='\t')

            if b == 'A':
                hdrs = m.columns

        except:
            print(f'batch {b} not present')
            m = pd.DataFrame(columns=hdrs)

            m['well'] = gt.well_range('A01', 'H12')
            m['type'] = 'empty'

        m['well'] = m['well'].copy().apply(lambda w: calc384well(w, b))
        m['batch'] = b

        fullplate = pd.concat([fullplate, m], axis=0)

        namelist.append(os.path.split(file)[-1].split('_')[0])

    fullplate.loc[fullplate['type']=='empty', 'batch'] = np.nan

    outdir = os.path.split(flist[0])[0]

    fullplate.sort_values('well', inplace=True)

    name = namelist[0]

    if len(set(namelist)) != 1:
        print("base plate names didn't agree")
        name += '_merged'

    fullplate.to_excel(os.path.join(outdir, name + '.xlsx'), index=False)


def batch_summary(file):
    """ summarizes identities of each batch in a plate map, one well per batch  """
    m = gct.openmap(file)
    batches = m['batch'].dropna().unique()
    res = []
    for b in batches:
        res.append(gt.hsub(m, {'batch': b}).iloc[3])
    res = pd.concat(res, axis=1)
    return res.T


def convert_to_txt(path):
    fl = gt.get_flist(path, '.xlsx')
    opath = '/Users/WRB/Desktop/newmaps/'
    for f in fl:
        m = pd.read_excel(f)
        shn = gt.get_shn(f).strip('.xlsx')
        outpath = os.path.join(opath, shn +'.txt')
        m.to_csv(outpath, sep='\t', index=False)


def list_fails(path):
    fl = gt.get_flist(path)
    for f in fl:
        if f.endswith('.txt'):
            m = pd.read_table(f, index_col=False)
        elif f.endswith('.xlsx'):
            m = pd.read_excel(f)
        fw = m[m['type']=='failed']['well'].values
        shn = gt.get_shn(f)
        print(shn, fw)


def compare_files(path1, path2, show=False):
    """ only works for text files, compares to see if identical """
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
    """ from two folders, compare to see whether all identically named maps are identical """
    fl1, fl2 = gt.get_flist(path1), gt.get_flist(path2)
    for file in fl1:
        shn = gt.get_shn(file)
        other = os.path.join(path2, shn)
        if os.path.isfile(other):
            compare_files(file, other, show=show)


def ctup(df, arg_dict, col=':', u=True, lst=False):
    """ count up, check to see how many entries in df satisfy the argdict """
    if lst is True and u is True:
        r = gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna().unique()
    elif lst is False and u is True:
        r = len(gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna().unique())
    elif lst is True and u is False:
        r = gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna()
    elif lst is False and u is False:
        r = len(gt.hsub(df, arg_dict=arg_dict).loc[:,col].dropna())
    return r


def plate_map_vis(myseries, cmap='dflt', path='dflt'):
    """ just translate directly into a dict or overwrite to use?
    returns just the array to plot """
    num_cats = len(myseries.unique())
    cat2num = dict(zip(myseries.unique(), range(num_cats)))
    data = myseries.apply(lambda x: cat2num[x])
    if cmap == 'dflt':
        if num_cats < 10:
            cmap = 'tab10'
            maxcats = 10
        else:
            camp = 'tab20'
            maxcats = 20
    if path is 'dflt':
        outpath = gt.dflt_outpath()
    else:
        outpath = path
    plottools.plot_plateplot(data, outpath=outpath, name=myseries.name, label=data.name, ncats=maxcats, cmap=cmap, clrbar=cat2num)


def summarize_doses(h):
    # returns series of names and lists of their unique doses as values
    h.sort_values(['batch', 'name','dose'])
    try:
        res = gt.hsub(h, {'type':'test'}).sort_values('dose')['dose'].groupby(h['name']).unique()
    except:
        res = gt.hsub(h, {'type': 'test'})['dilution'].groupby(h['name']).unique()
    return res



def split_tabs(file):
    """ take excel file with multiple maps in each tab and split out into
    individual excel files per tab """
    outpath = os.path.split(file)[0]
    pname = gt.get_shn(file).split('.')[0]
    xl = pd.ExcelFile(file)
    sheets = xl.sheet_names
    for sh in sheets:
        m = xl.parse(sh)
        # extras specific to this set
        m['batch'] = 'A'
        sh_name = sh.replace('_','-')
        savepath = os.path.join(outpath, pname + '-' + sh_name)
        m.to_excel(savepath + '.xlsx', index=False)



def check_maps(path, compare=True, img=True, v=True, filt=True):
    """ looks through .txt and .xlsx maps in a directory and summarizes their content and relationship with each other,
    as well as generating plate visualizations of type and batch for each plate. V for verbose, lists names and doses per plate
     if filter is true, just observe 6character name excel files, otherwise consider all files"""

    outpath = os.path.join(path, 'map summary')
    try:
        os.mkdir(outpath)
    except:
        pass

    plot_fields = ['type', 'batch']
    #plot_fields = ['type']
    # add checks to add batch and dose if present

    pert_dict, map_list = {}, {}
    wellpert_dict = {}

    flist = gt.get_flist(path, ext='.xlsx')
    if filt is True:
        flist = [x for x in flist if len(os.path.split(x)[-1]) == 11]
    if len(flist) == 0:
        flist = gt.get_flist(path, ext='.txt')
        flist = [x for x in flist if len(os.path.split(x)[-1])==10]

    if v is True:
        print('flist = ', flist)

    awells = gt.get_awells()
    composition = pd.DataFrame(columns=['wells #', 'test #', 'doses', 'dose/trt', '# names','vehicle #','poscon #','poscons'])

    for f in flist:
        pname = gt.get_shn(f).split('.')[0]
        print(pname)
        if f.endswith('.xlsx'):
            m = pd.read_excel(f)
        elif f.endswith('.txt'):
            m = pd.read_table(f, index_col=False)
        m.sort_index(inplace=True)
        batches = m['batch'].dropna().unique()

        if any([('dose' in x) or ('dilution' in x) for x in m.columns]):
            dose_field = [x for x in m.columns if (('dose' in x) or ('dilution' in x))][0]
        else:
            dose_field = None

        headers = {'wells #': lambda x: len(x.index),
                   'test #': lambda x: ctup(x, {'type': 'test'}, 'well')}
        if dose_field is not None:
            headers.update({
                'doses': lambda x: ctup(x, {'type': 'test'}, dose_field),
                'dose/trt': lambda x: gt.hsub(m, {'type': 'test'})[dose_field].groupby(m['name']).unique().apply(
                    lambda x: len(x)).mean()})
        elif dose_field is None:
            headers.update({'doses': 'na', 'dose/trt': 'na'})
        headers.update({
            '# names': lambda x: ctup(x, {'type': 'test'}, 'name'),
            'vehicle #': lambda x: ctup(x, {'type': 'vehicle'}, 'well'),
            'poscon #': lambda x: ctup(x, {'type': 'poscon'}, 'well'),
            'poscons': lambda x: ctup(x, {'type': 'poscon'}, 'name', lst=True)})

        summary = pd.DataFrame(columns=headers)

        # check wells for full plate
        well_result = set(awells) - set(m['well'].values)

        if len(well_result) != 0:
            print('{} wells error, {} entries'.format(pname, len(m.index)))

        if v is True:
            print(gt.hsub(m, {'type':'test'})['name'].dropna().unique())
            try:
                doselist = gt.hsub(m, {'type':'test'})[dose_field].dropna().unique()
                print(doselist)
            except:
                print('no dose column')
                pass

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
                except (KeyError, TypeError):
                    summary.loc[entry, k] = 'na'

        composition = pd.concat([composition, summary])

        if img is True:
            for pf in plot_fields:
                plot_series = m[pf]
                if len(plot_series.dropna().unique()) > 1:
                    plot_series.name = pname + ' ' + pf
                    plate_map_vis(plot_series, path=outpath)

    composition.to_excel(os.path.join(outpath, 'batch_composition.xlsx'))

    if compare is True:
        same_plates = gt.overlap_matrix(wellpert_dict.values(), wellpert_dict.keys())
        name_overlap = gt.overlap_matrix(pert_dict.values(), pert_dict.keys())
        name_overlap.to_excel(os.path.join(outpath, 'name_overlaps.xlsx'))
        same_plates.to_excel(os.path.join(outpath, 'well-name_overlaps.xlsx'))





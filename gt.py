#!/usr/bin/python
""" this is a general tools script which contains many common handy functions """

import os, csv, glob, string, shutil, gt, gct
import skyline
import pandas as pd
import numpy as np
import collections as cll
import matplotlib.pyplot as plt
import seaborn as sns
try:
    import pyperclip
except:
    pass


def get_well_reps(h, well, cats, df=False):
    """ return list of well addresses of other wells in the header file matching the
    provided categories as the passed well. n=name, d=dose, c=cell, b=batch
    returns list of addrs unless df is True, then it passes bach header dataframe of those wells"""
    args = []
    if len(well) == 3:
        well = h.index[0][:-3] + well
    if 'n' in cats:
        args.append('name')
    if 'd' in cats:
        args.append('dose')
    if 'b' in cats:
        args.append('batch')
    if 'c' in cats:
        args.append('cell')

    argdict = {}
    try:
        mywell = h.loc[well]
    except KeyError:
        print(f'{well} well not found in index')
        return 'empty'
    print(mywell)
    for a in args:
        argdict[a] = mywell[a]

    matches = gt.hsub(h, argdict)
    if df is True:
        return matches
    else:
        return list(matches.index.values)


def splitpaths(pathlist, ext):
    pathlist = pathlist.replace(ext, ext+'å').replace('Macintosh HD', '')[:-1]
    pathlist = pathlist.split('å')
    return pathlist


def get_genes(genes, df=None):
    """ returns genes, with shortcuts of 'test1', 'test2' which is 12, 'test100' will plot first 100
     genes, and 'all' will return all genes in index """
    # check for gene arguments, and fill 'test' and 'all'
    if genes is 'test1':
        genes = ['200678_x_at']
    elif genes is 'test2':
        genes = ['121_at', '219888_at', '218245_at', '206501_x_at', '203154_s_at',
                 '201614_s_at', '209682_at', '202324_s_at', '209603_at',
                 '200060_s_at', '202123_s_at', '201579_at']
    elif genes is 'test100':
        genes = df.index.values[50:150]
    elif genes is 'all':
        genes = df.index.values
    return genes


def consensus(ds, name='dflt'):
    """ merges ZScore instances in a passed dataframe into a conservitive consensus - min abs value
    as long as all reps agree on direction, otherwise zero. uses a fancy np.select method to set values
    if name is 'dflt' the column name will be: FPA001:K03-N04_(3)
    otherwise if name is 'first' name will just be first well FPA001:K03 (better for header integration) """
    choices = [ds.min(axis=1), ds.max(axis=1)]
    conds = [(ds > 0).all(axis=1), (ds < 0).all(axis=1)]
    newdf = pd.Series(data=np.select(conds, choices, default=0), index=ds.index)
    if name is 'dflt':
        newdf.name = ds.columns[0] + '-' + ds.columns[-1].split(':')[-1] + '_(' + str(len(ds.columns)) + ')'
    elif name is 'first':
        newdf.name = ds.columns[0]
    return newdf


def create_desc_label(h, order='dflt',):
    if order is 'dflt':
        try:
            i = h['dose']
            h['desc'] = h['batch'] + '-' + h['name'] + '-' + h['dose']
            # order = ['batch', 'name', 'dose']
        except KeyError:
            h['desc'] = h['batch'] + '-' + h['name']
            # order = ['batch', 'name']
    else:
        h['desc'] = h[order].apply(lambda x: '-'.join(x), axis=1)
    return h['desc']


def copyout(path, term, outpath='dflt', rc=True):
    """ local copy files function, specify destination path, searchterm with wilcards, oupath and recursive
    aruments assumed """
    fl = globit(path, term, rc=rc)
    if outpath is 'dflt':
        outpath = dflt_outpath()
    for f in fl:
        shutil.copy(f, outpath)


def globit(path, term, rc=True):
    """ a wrapped recursive glob function """
    fl = glob.glob(path + '/**/*' + term, recursive=True)
    return fl


def stack_files(path, orient='vert'):
    dlist = []
    fl = get_flist(path, ext='.txt')
    for f in fl:
        dat = pd.read_table(f, index_col=0)
        dlist.append(dat)
    if orient is 'vert':
        myaxis = 0
    else:
        myaxis = 1
    data = pd.concat(dlist, axis=myaxis)
    fn = os.path.basename(os.path.dirname(path)) + '_joined.txt'
    outpath = os.path.join(path, fn)
    data.to_csv(outpath, sep='\t')


def gather_rows(path, searchstring, ext='all', save=False):
    """" extract rows from files in directory matching string, with optional file extension """
    flist = get_flist(path, ext)
    results = []
    if isinstance(searchstring, str):
        searchstring = [searchstring]
    for file in flist:
        with open(file, 'r') as f:
            for line in f:
                if any([x in line for x in searchstring]):
                    results.append(os.path.split(file)[-1] + '\t' + line)

    if save is True:
        outpath = os.path.join(path, searchstring[0] + '_rows.txt')
        with open(outpath, 'w', newline='') as outf:
            for line in results:
                outf.write(line)
    else:
        return results


def dflt_outpath(fldr_name='foo', path='dflt', fn=None):
    if path is 'dflt':
        path = gt.check_desktop()
    path = path + fldr_name
    try:
        os.makedirs(path)
    except OSError:
        pass
    if fn:
        path = os.path.join(path, fn)
    return path


def overlap_matrix(mysets, labels):
    """ pass in a list of lists/sets of identifiers, and the corresponding names of those
    collections, and get in return a symmetric dataframe with pairwise count of overlaps amongst groups"""
    results = pd.DataFrame(index=labels, columns=labels)
    for cohort, label in zip(mysets,labels):
        for comp, clabel in zip(mysets,labels):
            res = len(set(cohort) & set(comp))
            results.loc[label, clabel] = res
    return results


def isnum(s):
    """ quick custom test if an object is a number, aka int function can work on it """
    try:
        int(s)
        return True
    except ValueError:
        return False
    except TypeError:
        return False


def dosub(d, h, arg_dict, name='dflt'):
    """ same as dsub, except only returns the subset dataframe without the accompanying subset header """
    hs = hsub(h, arg_dict)
    ds = d[hs.index.values].copy()
    if name == 'dflt':
        try:
            ds.name = d.name + '_sub'
        except AttributeError:
            pass
    else:
        ds.name = name
    return ds


def dsub(d,h, arg_dict, name='dflt'):
    """ pass the dataframe and header to return the subset header as well as corresponding df slice. retrns ds, hs """
    hs = hsub(h, arg_dict)
    ds = d[hs.index.values].copy()
    if name == 'dflt':
        try:
            ds.name = d.name + '_sub'
        except AttributeError:
            pass
    else:
        ds.name = name
    return ds, hs


def hsub(h, arg_dict, sep=False):
    """ takes in dictionary of {param: value} to scrape folder and return data meeting criteria,
    dictionary value may be a list, in which case any value of list is acceptable.
    returns the filtered header file. if d=dataframe is specified returns ds, hs as pre-filtered data too
    if sep is true, then a passed list will result in a list of the results rather than
    a combined result of all of them, so a list of lists"""
    sh = h
    for c, v in arg_dict.items():
        if isinstance(v, str):
            sh = sh[sh[str(c)] == str(v)].copy()
        else:
            try:
                if sep is True:
                    sh = [sh[sh[str(c)] == str(x)] for x in v].copy()
                elif sep is False:
                    sh = sh[sh[str(c)].isin([str(x) for x in v])].copy()
            except:
                print('error with filter dictionary value')
    if len(sh.index.values) == 0:
        print('no wells in selection', arg_dict)
    return sh


def breakdown(df,h,cats, dic=True, genes=None):
    """ takes a dataframe and header and the categories to break down by 'b' batch, 'c' cell, 'n' name, 'd' dose.
    returns a dictionary with the key as the description and the dataframe as the value.
    'w' is also supported as breakdown by well - useful for many plates with identical layout

    if dic is True a dictionary is returned, with a key title and dataframe value
    if dic is False then list is returned, of tuples with dataframe and header

    """

    if genes is not None:
        genes = get_genes(genes)
        df = df.loc[genes]

    if 'd' in cats:
        try:
            dose_col = [x for x in h.columns if 'dose' in x or 'dilution' in x][0]
        except IndexError:
            print('dose column error')
    else:
        dose_col = None

    vd = cll.OrderedDict()
    subs = []

    cd = {'c': 'cell',
          'b': 'batch',
          'd': dose_col,
          'n': 'name',
          'w': 'well',
          'p': 'plate'}

    clist = []

    for c in cats:
        try:
            clist.append(cd[c])
        except IndexError:
            print('error, more than 3 categories')

    cat1 = clist[0]
    group1 = sorted(h[cat1].dropna().unique())
    for e1 in group1:
        argdict = {cat1: e1}
        try:
            cat2 = clist[1]
            for e2 in sorted(gt.hsub(h, {cat1: e1})[cat2].dropna().unique()):
                argdict.update({cat2: e2})
                try:
                    cat3 = clist[2]
                    for e3 in sorted(gt.hsub(h, {cat1: e1, cat2: e2})[cat3].dropna().unique()):
                        argdict.update({cat3: e3})
                        hdr = f'{e1}-{e2}-{e3}'
                        if dic is True:
                            vd.update({hdr: gt.dosub(df,h, argdict, name=hdr)})
                        else:
                            subs.append(gt.dsub(df,h, argdict, name=hdr))
                except IndexError:
                    hdr = f'{e1}-{e2}'
                    if dic is True:
                        vd.update({hdr: gt.dosub(df, h, argdict, name=hdr)})
                    else:
                        subs.append(gt.dsub(df, h, argdict, name=hdr))
        except IndexError:
            hdr = f'{e1}'
            if dic is True:
                vd.update({hdr: gt.dosub(df, h, argdict, name=hdr)})
            else:
                subs.append(gt.dsub(df, h, argdict, name=hdr))

    if dic is True:
        return vd
    else:
        return subs


def assemble_consensus(df, h, cats, ccs=True, plot=False, skyl=False, n=None, save=False, test=False):
    """ tool to assemble replicate zscore consensus, pass df, header and the breakdown categories 'nd' for instance
    will return the consolidated df and header file

    ccs will calculate the zscore correlation of replicates, and insert that into header df
    plot will use seaborn pairplot to visualize the calculated rep correlations above
    skyl controls skyline plot generation, can be True to plot all ind reps plus consensus
    n argument is a limiter to only consider treatments with enough replicates, including into consensus gct!!
    save will save the consensus gct file
    """

    if isinstance(df, str):
        df, h = gct.extractgct(df)

    outpath = gt.dflt_outpath(fldr_name='output figs')
    pname = df.name
    try:
        os.mkdir(os.path.join(outpath, pname))
    except:
        pass

    outpath = os.path.join(outpath, pname)

    subs = breakdown(df, h, cats, dic=False)

    con_data = pd.DataFrame(index=df.index)
    if ccs is True:
        con_header = pd.DataFrame(index=np.concatenate([h.columns.values,['corr', 'all ccs']]))
    else:
        con_header = pd.DataFrame(index=h.columns)

    for ds, hs in subs:
        if n is not None:
            if len(ds.columns) < n:
                print('not enough reps', hs.iloc[0])
                continue

        c = consensus(ds, name='first')
        con_data = pd.concat([con_data, c], axis=1)

        new_annot = hs.iloc[0,:].copy().T
        new_annot.well = hs['well'].values
        new_annot.addr = hs['addr'].values

        if ccs is True:
            corrs = []
            for i in range(len(ds.columns)):
                for j in range(1 + i, len(ds.columns)):
                    corrs.append(round(ds.iloc[:, i].corr(ds.iloc[:, j], method='pearson'), 2))
            if len(corrs) == 0:
                print('corrs = na')
                print(hs.iloc[0].values)
                new_annot['corr'] = 'na'
                new_annot['all ccs'] = 'na'
            elif len(corrs) == 1:
                new_annot['corr'] = corrs
                new_annot['all ccs'] = corrs
            else:
                new_annot['corr'] = round(np.percentile(corrs, 75), 2)
                new_annot['all ccs'] = corrs

        if plot is True:
            ds.columns = [x + ' - ' + hs.loc[x]['batch'] for x in ds.columns]
            ax = sns.pairplot(ds)
            myoutpath = os.path.join(outpath, 'rep zs scatter')
            try:
                os.mkdir(myoutpath)
            except:
                pass
            plt.savefig(os.path.join(myoutpath, h.plate[0] + '-' + ds.name + '.png'))
            plt.close()

        con_header = pd.concat([con_header, new_annot], axis=1)

        if skyl is True:
            myoutpath = os.path.join(outpath, 'skyline')
            try:
                os.mkdir(myoutpath)
            except:
                pass
            name = hs.iloc[0]['name'] + '-' + hs.iloc[0]['dose']
            name = name.replace('.', ',')
            title = pname + '-' + name
            myoutpath = os.path.join(myoutpath, title)
            skyline.new_skyline(ds, title=title, outpath=myoutpath)

        if test is True:
            break

    con_header = con_header.T

    if save is False:
        return con_data, con_header
    elif save is True:
        gct.save_headergct(con_data, gt.dflt_outpath(fn=df.name+'_consensus.gct'), con_header)


def txt2list(file):
    """ assemble a list from text file with variable items per line, sep='\t'
        returns list with sublists of variable length. If only one item per line
        then a single list with each line as one string is delivered. Otherwise list of lists """
    mylist = []
    with open(file, 'rU') as f:
        for line in f:
            n = line.split('\t')
            mylist.append([x.strip() for x in n])
    if max([len(x) for x in mylist]) == 1:
        mylist = [x[0] for x in mylist]
    return mylist


def get_shn(file):
    """ return the plate name (first component of filename split by '_') """
    shn = os.path.split(file)[-1].split('_', 1)[0]
    return shn


def get_flist(path, ext='all', shn=False):
    """ assemble list of absolute paths of files in top level of folder
            with specified file extension, if shn=true return dictionary of
            with key as shortname and value as full file path in addition to flist"""
    f_list = []
    try:
        for input_file in os.listdir(path):
            fullfilename = os.path.abspath(os.path.join(path, input_file))
            # exclude hidden files and filter for extensions
            if input_file[0] != '.' and input_file[0] != '$':
                if ext is 'all':
                    f_list.append(fullfilename)
                elif ext is not 'all':
                    if input_file.endswith(ext):
                        f_list.append(fullfilename)
        # build out {shn: fullname} dictionary
        if shn is True:
            shnd = {}
            for f in f_list:
                shnd[get_shn(f)] = f
            return f_list, shnd
        else:
            return f_list
    except TypeError:
        print('error, file extension for file list builder not provided')


def check_desktop():
    """ grab the current desktop working directory for whichever machine in use,
            need to add additional options to the list if desired """
    dtops = ['/Volumes/WRBHDD/wrb/Desktop/', '/Users/wrb/Desktop/', '/home/wrb/Desktop']
    for d in dtops:
        if os.path.exists(d):
            return d


def tolist(mystring, spl='_', uniq=False):
    """ returns list of plate short names from absolute paths
        default list, if uniq=True is passed then a unique set is returned"""
    # convert terminal escaped spaces into underscore
    try:
        lst = mystring.replace('\ ', '_')
    except AttributeError:
        lst = mystring
    lst = lst.split()
    lst = [x.strip() for x in lst]
    lst = [os.path.split(x)[1] for x in lst]
    if spl is not None:
        lst = [x.split(spl)[0] for x in lst]
    if uniq is False:
        return lst
    elif uniq is True:
        return set(lst)


def addr_id_df(df, p=None):
    """ convert dataframe into plate:well format column headers, if l is true just convert list"""
    if p is not None:
        if ':' not in df[0]:
            df = [p + ':' + w for w in df]
        else:
            print('already composite addr')
        return df
    if len(df.columns[1]) == 3:
        df.columns = df.shortname.split('-')[0] + ':' + df.columns
    elif len(df.columns[1]) != 3:
        if '_' in df.columns[4]:
            pname = df.columns.str.split('_', 1).str[0]
            wid = df.columns.str.split(':', 1).str[1]
            df.columns = pname + ':' + wid
        else:
            pass
    return df


def savelist(l, path):
    """ save the list to text file in given path, one item per line """
    with open(path, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        for el in l:
            csvwriter.writerow([el])


def clip(l):
    """ using pyperclip, copy a passed list to clipboard, one per line """
    s = str(l)
    chars = ['[', ']', '{', '}', ',', "'", '"']
    for c in chars:
        s = s.replace(c, '')
    s = s.replace(' ', '\n')
    s = s.replace('\n\n', '\n')
    if pyperclip:
        pyperclip.copy(s)
        print('list copied to clipboard (one element per line)')
    else:
        print('sorry, pyperclip not installed')


def tally_failed_wells(path):
    """ searches through folder and tallys count of failed wells in each plate """
    fl = glob.glob(path + '/**/*badwell_id*', recursive=True)
    wc = cll.Counter()
    for f in fl:
        with open(f, 'r') as infile:
            wells = infile.readline().split(',')
            wells = [x.strip() for x in wells]
            wc.update(wells)
    # pickle.dump(wc, open('wc.p', 'wb'))
    wc = sorted(wc.items(), key=lambda x: x[1], reverse=True)
    for item in wc:
        print(item)


def count_mode_failures(path):
    """ counts how many wells fail for indicated failure mode only
    # and how many fail with that and one other, versus the total """
    fl = glob.glob(path + '/**/*QC_fail*', recursive=True)
    # the header of the failuremode in txt file
    fmode = 8
    fails = {}
    for file in fl:
        with open(file, 'rU') as f:
            plate = os.path.split(file)[1].split('_')[0]
            linereader = csv.reader(f, delimiter='\t')
            fails[plate] = [0, 0, 0, 0]
            next(linereader)
            line = next(linereader)
            col = line.index('failedQC')
            next(linereader)
            next(linereader)
            for line in linereader:
                qcf = [int(x.strip()) for x in line[col].split()]
                # count if failure mode occurs for each well
                if fmode in qcf and len(qcf) == 1:
                    fails[plate][0] += 1
                elif fmode in qcf and len(qcf) == 2:
                    fails[plate][1] += 1
                else:
                    fails[plate][2] += 1
                # increment overall
                fails[plate][3] += 1
    # combine results and save
    fp = pd.DataFrame(fails)
    fp = fp.T
    op = os.path.join(path, 'fail_summary.txt')
    fp.to_csv(op, sep='\t')


def separate_subset_folders(path, mylist, down=False, dest='dflt'):
    """ copy top level folders over from path to destination if the folders match any terms in
        mylist, all folders transferred over as-is. if not found, printed """

    if dest is 'dflt':
        dest = gt.dflt_outpath()

    for st in mylist:
        if down is False:
            fl = glob.glob(path + st + '*')
        elif down is True:
            fl = glob.glob(path + '*/' + st)

        cplist = [x for x in fl if os.path.isdir(x)]

        try:
            dirpath = cplist[0]
            bn = os.path.basename(dirpath)
            try:
                shutil.copytree(dirpath, os.path.join(dest, bn))
            except FileExistsError:
                pass
        except IndexError:
            print(st, ' not found')


def separate_subset(path, dest, mylist, st1=None, st2=None, d=False):
    """ copy files from path folder to dest folder, matching items in list
            st1 and st2 are search terms (in order), d option deletes files
            only works to copy over flat files, will not preserve directory struct"""
    new, cplist = [], []
    if st1 is not None:
        new.append(st1)
    if st2 is not None:
        new.append(st2)
    new = '*'.join(new)
    slist = path + '**/*' + new
    if slist.endswith('*'):
        pass
    else:
        slist = slist + '*'
    print(slist)

    fl = glob.glob(slist, recursive=True)

    for f in fl:
        if any([x in f for x in mylist]):
            cplist.append(f)
    for f in cplist:
        fn = os.path.split(f)[1]
        print(f)
        fdest = os.path.join(dest, fn)
        shutil.copy(f, fdest)

    if d is True:
        for f in cplist:
            os.remove(f)


def get_awells(pos=True, ref=True, proc=True, empt=True):
    """ returns list of 384 three character well IDs """
    awells = []
    rows = string.ascii_uppercase[0:16]
    cols = range(1, 25)
    for l in rows:
        for n in cols:
            # join lettrs and nums with 2 characater num format
            awells.append(l + str('{:02d}'.format(n)))
    if proc is False:
        pos, ref, empt = False, False, False
    if pos is False:
        awells.remove('B01')
    if ref is False:
        awells.remove('B02')
        awells.remove('A02')
    if empt is False:
        awells.remove('A01')
    return awells


def wells_range(well_list):
    """ returns list of 3char ids of all well range tuples in the list, can string multiple together
    well_list = [('A23','E23'),('B01','F06')]"""
    wells = []
    if 'top' in well_list or 'bot' in well_list:
        if 'top' in well_list:
            wells.extend(well_range('C11', 'F14'))
        if 'bot' in well_list:
            wells.extend(well_range('K11', 'N14'))
        return wells
    else :
        for w in well_list:
            wells.extend(well_range(w[0], w[1]))
        print('length: ', len(wells))
        return wells


def well_range(upper_left, lower_right):
    """ returns list of 3char ids of all wells within provided rectangle coords """
    well_list = []
    alpha = string.ascii_uppercase
    startlet, startcol = upper_left[0], int(upper_left[1:])
    stoplet, endcol = lower_right[0], int(lower_right[1:])
    startpos = alpha.index(startlet)
    stoppos = alpha.index(stoplet)
    for l in alpha[startpos:stoppos + 1]:
        for i in range(startcol, endcol + 1):
            well = l + '{:02}'.format(i)
            well_list.append(well)
    return well_list


def main():
    print('no main')


if __name__ == '__main__':
    main()

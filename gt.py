#!/usr/bin/python
""" this is a general tools script which contains many common handy functions """

import os, csv, glob, string, shutil, gt
import pandas as pd
import numpy as np
import collections as cll
try:
    import pyperclip
except:
    pass

def get_genes(genes, df=None):
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


def consensus(ds):
    new = pd.DataFrame(index=ds.index)
    choices = [ds.min(axis=1), ds.max(axis=1)]
    conds = [(ds > 0).all(axis=1), (ds < 0).all(axis=1)]
    new['consensus']= np.select(conds, choices, default=0)
    return new


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


def dosub(d, h, arg_dict):
    """ same as dsub, except only returns the subset dataframe without the accompanying subset header """
    hs = hsub(h, arg_dict)
    ds = d[hs.index.values].copy()
    try:
        ds.name = d.name + '_sub'
    except AttributeError:
        pass
    return ds


def dsub(d,h, arg_dict):
    """ pass the dataframe and header to return the subset header as well as corresponding df slice. retrns ds, hs """
    hs = hsub(h, arg_dict)
    #print(hs.index[:5])
    ds = d[hs.index.values].copy()
    try:
        ds.name = d.name + '_sub'
    except AttributeError:
        pass
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


def breakdown(df,h,cats):
    """ takes a dataframe and header and the categories to break down by 'b' batch, 'n' name and 'd' dose.
    will return a list of dataframe/header pairs for each respective unique combo of the given batches.
    useful to subset a df and then plot the remaining df along one of these categories """
    cd = {}

    if 'd' in cats:
        try:
            dose_col = [x for x in h.columns if 'dose' in x][0]
        except:
            print('dose column error')

    subs = []


    if ''.join(sorted(cats)) is 'bdn':
        for b in h.batch.dropna().unique():
            for n in h.name.dropna().unique():
                for d in h.dose_col.dropna().unique():
                    subs.append(dsub(df,h,{'batch':b, 'name':n, dose_col:d}))
    elif ''.join(sorted(cats)) is 'bn':
        for b in h.batch.dropna().unique():
            for n in h.name.dropna().unique():
                subs.append(dsub(df,h,{'batch':b, 'name':n}))
    elif cats is 'b':
        for b in h.batch.dropna().unique():
            subs.append(dsub(df,h,{'batch':b}))
    elif cats is 'n':
        for n in h.name.dropna().unique():
            subs.append(dsub(df,h,{'name':n}))
    elif cats is 'd':
        for d in h.dose_col.dropna().unique():
            subs.append(dsub(df,h,{dose_col:d}))

    return subs


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
    dtops = ['/Volumes/WRBHDD/wrb/Desktop/', '/Users/wrb/Desktop/']
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

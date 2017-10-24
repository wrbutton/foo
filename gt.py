#!/usr/bin/python
""" this is a general tools script which contains many common handy functions """

import os, csv, glob, string, shutil
import pandas as pd
import collections as cll
try:
    import pyperclip
except:
    pass


def hsub(h, arg_dict):
    """ takes in dictionary of {param: value} to scrape folder and return data meeting criteria,
    dictionary value may be a list, in which case any value of list is acceptable.
    returns the filtered header file"""
    sh = h
    for c, v in arg_dict.items():
        if isinstance(v, str):
            sh = sh[sh[c] == v]
        elif isinstance(v, list):
            sh = sh[sh[c].isin(v)]
        else:
            print('error with filter dictionary value')
    if len(sh.index.values) == 0:
        print('no wells in selection', arg_dict)
    return sh


def txt2list(file):
    """ assemble a list from text file with variable items per line, sep='\t'
        returns list with sublists of variable length. If only one item per line
        then a single list with each line as one string is delivered. Otherwise list of lists """
    mylist = []
    with open(file, 'rU') as f:
        for line in f:
            n = line.split('\t')
            mylist.append([x.strip() for x in n])
    if max([len(x.split('\t')) for x in mylist]) == 1:
        mylist = [x[0] for x in mylist]
    return mylist


def get_shn(file):
    """ return the plate name (first component of filename split by '_') """
    shn = os.path.split(file)[-1].split('_', 1)[0]
    return shn


def get_flist(path, ext, shn=False):
    """ assemble list of absolute paths of files in top level of folder
            with specified file extension, if shn=true return dictionary of
            with key as shortname and value as full file path in addition to flist"""
    f_list = []
    try:
        for input_file in os.listdir(path):
            fullfilename = os.path.abspath(os.path.join(path, input_file))
            # exclude hidden files and filter for extensions
            if input_file[0] != '.':
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


def addr_id_df(df):
    """ convert dataframe into plate:well format column headers"""
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


def separate_subset(path, dest, mylist, st1=None, st2=None, d=False):
    """ copy files from path folder to dest folder, matching items in list
            st1 and st2 are search terms (in order), d option deletes files"""
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

    if d is True:
        for f in cplist:
            os.remove(f)
    else:
        for f in cplist:
            fn = os.path.split(f)[1]
            print(f)
            fdest = os.path.join(dest, fn)
            shutil.copy(f, fdest)


def get_awells():
    """ returns list of 384 three character well IDs """
    awells = []
    rows = string.ascii_uppercase[0:16]
    cols = range(1, 25)
    for l in rows:
        for n in cols:
            # join lettrs and nums with 2 characater num format
            awells.append(l + str('{:02d}'.format(n)))
    return awells


def well_range(startlet, stoplet, startcol, endcol):
    """ returns list of 3char ids of all wells within provided rectangle coords """
    well_list = []
    alpha = string.ascii_uppercase
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

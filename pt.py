#!/usr/bin/python

import gct, os, csv, glob, pickle
import shutil
import pandas as pd
import collections as cll


def txt2list(file):
    """ assemble a list from text file with variable items per line, sep='\t'
            returns list with sublists of variable length """
    mylist = []
    with open(file, 'rU') as f:
        for line in f:
            n = line.split('\t')
            mylist.append([x.strip() for x in n])
    return mylist


def get_shn(file):
    """ return the plate name (first component of filename split by '_') """
    shn = os.path.split(file)[-1].split('_',1)[0] 
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
        if shn == True:
            shnd = {}
            for f in f_list:
                shnd[get_shn(f)] = f
            return f_list, shnd
        else:
            return f_list
    except TypeError:
        print('error, file extension for file list builder not provided')


def getdeskpath():
    """ grab the current desktop working directory for whichever machine in use,
            need to add additional options to the list if desired """
    dtops = ['/Volumes/WRBHDD/wrb/Desktop/', '/Users/wrb/Desktop/']
    for d in dtops:
        if os.path.exists(d):
            return d


def tolist(string, s=False):
    """ returns list of plate short names from absolute paths
            default list, if s=True is passed then a set is returned"""
    # convert terminal escaped spaces into underscore
    l = string.replace('\ ', '_')
    l = l.split()
    l = [x.strip() for x in l]
    l = [os.path.split(x)[1] for x in l]
    l = [x.split('_')[0] for x in l]
    if s == False:
        return l
    elif s == True:
        return set(l)


def convert_to_addr(path):
    """ loop through gct files in folder and ensure that the column headers
            are in plate:well format, rewriting over files as necessary """
    flist = get_flist(path, '.gct')
    for f in flist:
        shn = get_shn(f).split('-')[0]
        print(shn)
        d, h = gct.extractgct(f)
        # convert df column headers to address
        d = addr_id_df(d)
        # if header file needs converting, do so
        if ':' not in h.index[4]:
            h.index = shn + ':' + h.index
        # strip down to basic plate name if legacy long format
        if '_' in h.index[2]:
            # use the dataframe.str method to apply regular commands
            pname = h.index.str.split('_', 1).str[0]
            wid = h.index.str.split(':', 1).str[1]
            h.index = pname + ':' + wid
        # save, overwriting file
        gct.save_headergct(d, f, h)


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


def platelist(string):
    d = getdeskpath()
    fp = os.path.join(d, 'platelist.txt')
    l = tolist(string)
    savelist(l, fp)


def savelist(l, path):
    with open(path, 'w') as out_file:
        csvwriter = csv.writer(out_file, delimiter='\t')
        for el in l:
            csvwriter.writerow([el])


def receipts(rcpts):
    l = rcpts.replace('\ ', '_').strip()
    l = rcpts.replace('\\', '')
    l = l.split()
    l = [x.strip() for x in l]
    l = [os.path.split(x)[1] for x in l]
    print(l)
    l = [float(x.split('_')[2]) for x in l]
    print(l)
    return sum(l)


def swap_index(df, old, new):
    df[old] = df.index.values
    df.set_index(new, inplace=True)
    df[new] = df.index.values
    return df


def edit_files(path):
    skiplines = 2
    flist = gct.get_flist(path, 'consensus.gct')
    for f in flist:
        shn = os.path.split(f)[1].split('_')[0]
        # define new file name here
        outpath = f.replace('consensus', 'consensus2')
        with open(f, 'rU') as in_file, open(outpath, 'w', newline='') as out_file:
            linereader = csv.reader(in_file, delimiter='\t')
            csvwriter = csv.writer(out_file, delimiter='\t')
            # skip down until reach place of interest
            for i in range(skiplines):
                csvwriter.writerow(next(linereader))
            # now modify lines of interest
            line = next(linereader)
            if ':' not in line[3]:
                line = [shn + ':' + x for x in line]
                line[0] = 'id'
            csvwriter.writerow(line)
            # then copy the rest overy
            for line in linereader:
                csvwriter.writerow(line)


def check_headers(path):
    skiplines = 2
    flist = gct.get_flist(path, '.gct')
    for f in flist:
        with open(f, 'rU') as in_file:
            linereader = csv.reader(in_file, delimiter='\t')
            for i in range(skiplines):
                next(linereader)
            line1 = len(next(linereader))
            line2 = len(next(linereader))
            for i in range(skiplines*4):
                next(linereader)
            line3 = len(next(linereader))
        result = line1 == line2 == line3
        print('ok!')
        if result == False:
            print(result, f)


def tally_failed_wells(path):
    fl = glob.glob(path + '/**/*badwell_id*', recursive=True)
    wc = cll.Counter()
    for f in fl:
        with open(f, 'r') as f:
            wells = f.readline().split(',')
            wells = [x.strip() for x in wells]
            wc.update(wells)
    #pickle.dump(wc, open('wc.p', 'wb'))
    wc = sorted(wc.items(), key=lambda x:x[1], reverse=True)
    for item in wc:
        print(item)


def count_mode_failures(path):
    # counts how many wells fail for indicated failure mode only
    # and how many fail with that and one other, versus the total
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


def comparewells(path, wids, mode='fc'):
    results = []
    flist = gct.get_flist(path, '.gct')
    for f in flist:
        d = gct.dfsubset(f, wids)
        if len(d.columns) < 2:
            continue
        if mode == 'fc':
            ddiff = d.iloc[:,1] / d.iloc[:,0]
        if mode == 'diff':
            ddiff = d.iloc[:,1] - d.iloc[:,0]
        ddiff.name = d.columns[0]
        ddiff.sort_values(ascending=False, inplace=True)
        results.append(ddiff)
    tc, bc = cll.Counter(), cll.Counter()
    for r in results:                                                 
        top = r.iloc[0:50].index.values
        tc.update(top)
        bot = r.iloc[-50:].index.values
        bc.update(bot)
    tc = [x for x in tc.items() if x[1] != 1]
    bc = [x for x in bc.items() if x[1] != 1]
    tc = sorted(tc, key=lambda x:x[1])
    bc = sorted(bc, key=lambda x:x[1])
    return results, bc, tc


def separate_subset(path, dest, mylist, st1=None, st2=None, d=False):
    ''' copy files from path folder to dest folder, matching items in list
            st1 and st2 are search terms (in order), d option deletes files'''
    new, cplist = [], []
    if st1 != None:
        new.append(st1)
    if st2 != None:
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

    if d == True:
        for f in cplist:
            os.remove(f)
    else:
        for f in cplist:
            fn = os.path.split(f)[1]
            print(f)
            fdest = os.path.join(dest, fn)
            shutil.copy(f, fdest)


def loud_wells(path):
    """ arg1 = path. looks through for gcts and makes counter for wells at the top"""
    flist = get_flist(path, '.gct')
    for f in flist:
        d, h = gct.extractgct(f)
        d = abs(d)
        p = d.sum(axis=0).sort_values(ascending=False)
        fn = get_shn(f)
        outpath = os.path.join(path, fn + 'promiscuity.txt')
        p.to_csv(outpath, sep='\t')

    wells = []
    numwells = 50
    promlist = get_flist(path, 'promiscuity.txt')
    for f in promlist:
        with open(f, 'rU') as file:
            for i in range(1,51):
                line = file.readline()
                well = line.split('\t')[0].split(':')[1]
                wells.append(well)
    
    cdct = cll.Counter(wells)

    cdct = pd.Series(cdct)
    cdct = cdct.sort_values(ascending=False)
    cdct.to_csv(os.path.join(path, 'loud_wells.txt'), sep='\t')


def get_awells():
  # returns list of 384 three character well IDs
  awells = []
  #rows = string.ascii_uppercase[0:16]
  rows = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"]
  cols = range(1,25)
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
    for l in alpha[startpos:stoppos+1]:
        for i in range(startcol, endcol+1):
            well = l + '{:02}'.format(i)
            well_list.append(well)
    return well_list

def comparetoplate(path, wids, mode='fc'):
    """ compare well ids to plate average, either by fold change or substraction
            input contains directory of gct files, wids in 3 char only, code will add
            the gct shortname to populate verbose well ids """
    results = []
    flist = gct.get_flist(path, '.gct')
    for f in flist:
        g = gct.Gct(f)
        d = g.build_dframe()
        av = d.mean(axis=1)
        l = [g.shortname + ':' + x for x in wids]
        print(l)
        d = d[d.columns[d.columns.isin(l)]]
        if len(d.columns) < 1:
            continue
        for w in l:
            if mode == 'fc':
                ddiff = d[w] / av
            if mode == 'diff':
                ddiff = d[w] - av
            ddiff.name = d.columns[0]
            ddiff.sort_values(ascending=False, inplace=True)
            results.append(ddiff)
    tc, bc = cll.Counter(), cll.Counter()
    for r in results:
        top = r.iloc[0:50].index.values
        tc.update(top)
        bot = r.iloc[-50:].index.values
        bc.update(bot)
    print(tc)
    tc = [x for x in tc.items() if x[1] != 1]
    bc = [x for x in bc.items() if x[1] != 1]
    tc = sorted(tc, key=lambda x:x[1])
    bc = sorted(bc, key=lambda x:x[1])
    return results, bc, tc


def main():
    print('no main')
if __name__ == '__main__':
    main()
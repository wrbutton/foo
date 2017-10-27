#!/usr/bin/python
'''
01/10/17 - added naming items in panel assembly with shortname
01/05/17 - added support to get names dictionary with well ids
12/11/16 - gct class object + helpers version 2.0

handles the interpretation of gct #1.3 files, the primary actions are
reading the header information and loading the data matrix into a 
pandas dataframe with wells as column ids and affy feature names as 
row indexes. 

will load header info if available, and store vehicle and poscon lists,
and place np.nan values over empty values. also stores all header lines
to be able to write along with modified data into a new gct-compatible 
object. 

supports a couple of additional seperate commands for assembling
a panel of multiple dataframes in a single object, and a quick 
method for taking the average or the median of a panel. 

'''
import pandas as pd
import numpy as np
import collections as cll
import os, csv, sys, pickle, pt

def gather_wells(path, wellfile):
    # reads from a text file ('MKD031:A03' one or many per line) to assemble a dataframe
    # from a collection of wells from a range of plates. dir points to gct directory
    fl = pt.get_flist(path, '.gct')
    data, h = pd.DataFrame(), pd.DataFrame()
    with open(wellfile, 'rU') as f:
        for line in f:
            plates = list(set([p.split(':')[0] for p in line.strip().split(',')]))
            wells = line.strip().split(',')
            # if only one plate, grab all wells from it
            if len(plates) <= 1:
                file = [s for s in fl if plates[0] in s]
                if len(file) > 1:
                    print('more than one plate match', file)
                d, sh = dfsubset(file[0], wells)
                d = d.sort_index()
                data = pd.concat([data, d], axis=1)
                h = pd.concat([h, sh], axis=0)
            # if multiple plates, just grab each well one at a time
            elif len(plates) > 1:
                for w in wells:
                    p = w.split(':')[0]
                    file = [s for s in fl if p in s]
                    if len(file) > 1:
                        print('more than one plate match', file)
                    d, sh = dfsubset(file[0], w)
                    d = d.sort_index()
                    data = pd.concat([data, d], axis=1)
                    h = pd.concat([h, sh], axis=0)
    return data, h


def dfsubset(path, wids):
    """ build a dataframe with only a selection of wells as samples, will create the dataframe in
    a memory friendly way, only loading those lines rather than opening the whole file to memory """
    g = Gct(path)
    g.get_headers()
    g.get_features()
    d, h = g.build_subset(wids)
    return d, h


def get_combined(path, ctype='median'):
    """ quick method to get median or mean of a panel """
    panel = load_panel(path, pan=True)
    if ctype == 'median':
        combined = panel.median(axis=0)
    elif ctype == 'mean':
        combined = panel.mean(axis=0)
    return combined


def load_panel(f_list, indv=False):
    # given folder path gets gct info and builds dframes for files
    # either merges into panel, or returns 2 dicts of gcts and dframes
    fnames, gd, dfd = [], {}, {}
    for i, file in enumerate(f_list):
        gn, dfn = 'g' + str(i), 'df' + str(i)
        gd[gn] = Gct(file)
        dfd[dfn] = gd[gn].build_dframe()
        fnames.append(gd[gn].shortname)
    if indv == True:
        return gd, dfd
    else:
        queue = zip(fnames, dfd.values())
        panel = pd.Panel({n: df for n, df in queue})
        return panel


def split_gct_batch(file):
    d, h = extractgct(file)
    batches = h['batch'].dropna().unique()
    if len(batches) < 2:
        print('no batches to split')
    for b in batches:
        hb = h[h['batch']==b]
        wells = hb.index.values
        db = d[wells]
        dirp = os.path.split(file)[0]
        fname = os.path.basename(file).split('_')[0]
        outpath = os.path.join(dirp, fname + b + '.gct')
        save_headergct(db, outpath, hb)


def get_collapsed(dflist, ctype='median'):
    """consolidate each dataframe in list and return consensus profile of each"""
    new_dfs = []
    for df in dflist:
        df.sort_index(inplace=True)
        title = df.shortname
        if ctype == 'median':
            profile = df.median(axis=1)
        elif ctype == 'mean':
            profile = df.mean(axis=1)
        profile.name = title
        new_dfs.append(profile)
    new_df = pd.concat(new_dfs, axis=1)
    return new_df


# incomplete, but to read a non-gct format matrix file into dframe
def save_simplegct(df, outpath):
    with open(outpath, 'w') as f:
        print('#1.3', file=f)
        print('{}\t{}\t0\t0'.format(len(df.index), len(df.columns)), file=f)
    df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    df.to_csv(outpath, mode='a', sep='\t', float_format='%.3f')


def save_headergct(df, outpath, h):
    with open(outpath, 'w') as f:
        wr = csv.writer(f, delimiter='\t')
        wr.writerow(['#1.3'])
        wr.writerow([len(df.index), len(df.columns), 0, len(h.columns)])
        wells = df.columns.values        
        row = ['well']
        row.extend(wells)
        wr.writerow(row)
        for col in h.columns:
            row = [col]
            vals = [h.ix[w, col] if w in h.index.values else 0 for w in wells]
            row.extend(vals)
            wr.writerow(row)
    with open(outpath, 'a') as f:
        df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        df.to_csv(f, sep='\t', header=False, float_format='%.3f')


def save_multiplateheadergct(df, outpath, h):
    # in process
    with open(outpath, 'w') as f:
        wr = csv.writer(f, delimiter='\t')
        wr.writerow(['#1.3'])
        wr.writerow([len(df.index), len(df.columns), 0, len(h.columns)])
        wells = df.columns.values
        row = ['id']
        row.extend(wells)
        wr.writerow(row)
        for col in h.columns:
            vals = []
            row = [col]
            for w in df.columns:
                try:
                    # try lookup by address and header column
                    val = h.loc[w, col]
                except KeyError:
                    print(w, 'keyerrer')
                    # look up by well instead of addr, only take first - RISKY!!
                    val = h[h['well']==w.split(':')[1]][col].values[0]
                vals.append(val)
            row.extend(vals)
            wr.writerow(row)
    with open(outpath, 'a') as f:
        df.apply(lambda x: pd.to_numeric(x, errors='ignore'))
        df.to_csv(f, sep='\t', header=False, float_format='%.3f')


# method to save a dafaframe object with gct file structure and headers
def save_dfgct(gct, df, outpath):
    with open(outpath, 'w') as f:
        for line in gct.hlines:
            print(line, end='', file=f)
    for i in range(gct.hcols):
        df.insert(0, 'foo' + str(i), 'nan')
    df.fillna('nan').to_csv(outpath, mode='a', sep='\t', float_format='%.2f')


def get_hd():
    hd = {'batch':'batch',     'det_well':'well',             'det_plate':'plate',
                'dose':'dose',         'concentration':'dose',    'pert_type':'type',
                'cell_id':'cell',    'pert_desc':'name',            'dilution':'dose', 'addr':'addr'}
    return hd


def extractgct(path):
    g = Gct(path)
    d = g.build_dframe()
    d.name = g.shortname
    try:
        h = g.fh
    except AttributeError:
        h = None
    return d, h


class Gct(object):
    # an object ossociated with extracting info from written gct files
    def __init__(self, infile):
        # pass the file argument as the filepath referenc
        self.file = infile
        self.shortname = os.path.split(infile)[-1].split('_',1)[0]

    def get_headers(self):
        # open file and extract the gct header information, storing dictionary values
        with open(self.file, 'rU') as in_file:
            line = in_file.readline()
            #check for .gct identity
            if line.strip() != '#1.3':
                sys.exit(self.file + 'not a .gct file')
            # add line to header lines list
            self.hlines = [line,]
            line = in_file.readline()
            self.hlines.append(line)
            # split and store gct dimensions from 2nd line
            hvals = [int(x) for x in line.split('\t')[:4]]
            # assign object variables for dimensions for easy retrieval
            self.drows = hvals[0]
            self.dcols = hvals[1]
            self.hcols = hvals[2]
            self.hrows = hvals[3]
            # !!! is this right? missing one row in some files
            self.skrows = hvals[3] + 2
            self.scol = hvals[2] + 1
            # grab well id's for matrix, if no header fields it's just line 3
            line = in_file.readline()
            self.addr = [x.strip() for x in line.split('\t')[self.scol:]]
            if self.hrows == 0:
                self.wells = [x.strip() for x in line.split('\t')[self.scol:]]
            # otherwise look through header lines and extract well row
            else:
                # for internal format gcts, internal metadata
                hd = get_hd()
                fheader = pd.read_csv(self.file, sep='\t', skiprows=2, nrows=self.hrows)
                fheader.set_index(fheader.columns[0], inplace=True)
                dropcols = [x for x in range(0,self.hcols)]
                fh = fheader.T
                fh.drop(fheader.columns[dropcols], inplace=True)
                # this is for the pipeline headers
                a = [x for x in fh.columns if x in hd.keys()]
                # this is for 'regular' headers
                b = [x for x in fh.columns if x in hd.values()]
                if len(a) > len(b):
                    fh = fh[a]
                    newnames = [hd[c] for c in fh.columns]
                    fh.columns = newnames
                    fh['addr'] = fh.index.values
                else:
                    fh = fh[b]
                fh.replace('na', np.NaN, inplace=True)
                fh.dropna(axis=1, how='all', inplace=True)
                fh.dropna(axis=0, how='all', inplace=True)
                #fh.str.replace({'type': {'trt_cp':'test', 'trt_poscon':'poscon',
                #                                             'ctl_vehicle':'vehicle', 'lma_x':'empty'}}, inplace=True)
                fh['type'].replace({'trt_cp':'test', 'trt_poscon':'poscon',
                                                             'ctl_vehicle':'vehicle', 'lma_x':'empty'}, inplace=True)
                self.fh = fh
                self.wells = fh.index.values
                self.vehicles = fh[fh['type'] == 'vehicle'].index.values
                self.poscons = fh[fh['type'] == 'poscon'].index.values
                self.empty = fh[fh['type'] == 'empty'].index.values
                self.test = fh[fh['type'] == 'test'].index.values
                try:
                    self.names = fh[fh['type'] == 'test']['name'].unique()
                except KeyError:
                    self.names = []
                self.batches = fh['batch'].unique()
                try:
                    self.doses = fh[fh.index.isin(self.test)]['dose'].unique()
                except KeyError:
                    pass

    def get_features(self):
        self.genes = []
        if not hasattr(self, 'skrows'):
            self.get_headers()
        with open(self.file, 'rU') as in_file:
            for i in range(self.skrows):
                in_file.readline()
            for line in in_file:
                self.genes.append(line.split('\t')[0].strip())

    def null_empty(self, data):
        # used to replace process control wells with 'nan' so they won't plot
        for w in self.empty:
            data.ix[:,w] = np.nan
        return data

    def build_dframe(self, nullempt=False):
        # open file and extract the gct header information, storing dictionary values
        self.get_headers()
        # assign the extracted array to data variable
        data = pd.read_csv(self.file, sep='\t', index_col=0, skiprows=self.skrows)
        data.drop(data.columns[range(0,self.hcols)], axis=1, inplace=True)
        data.columns = self.wells
        data.index.name = 'well'
        data.shortname = self.shortname
        # null out the process control wells
        if nullempt == True:
            data = self.null_empty(data)
        return data

    def build_sample_subset(self, wellids):
        # extract a subset of wells from a gct into a dataframe, without loading full
        # well ids should be supplied as full plate:well 9 character well ids
        self.get_headers()
        coli = []
        h = self.fh
        # try loading both the index and addr column to have a match for the column ids
        # for both single items as well as list for the col/row ids
        if isinstance(wellids, str):
            try:
                subh = h[h['addr']==wellids]
            except KeyError:
                subh = h[h.index==wellids]
        else:
            try:
                subh = h[h['addr'].isin(wellids)]
            except KeyError:
                subh = h[h.index.isin(wellids)]
        cols = subh.index.values
        for c in cols:
            coli.append(h.index.get_loc(c))
        coli = [x + (self.hcols+1) for x in coli]
        coli.insert(0, 0)
        data = pd.read_csv(self.file, sep='\t', index_col=0, skiprows=self.skrows+1,
                           usecols=coli, header=None)
        data.columns = cols
        print(data.columns)
        return data, subh

    def build_subset(self, ids):
        # extract a subset of wells from a gct into a dataframe, without loading full
        # well ids should be supplied as full plate:well 9 character well ids
        self.get_headers()
        h = self.fh
        cols, coli = [], []
        rows = []
        # wrap a single id in a list so can iterate with the same code
        if isinstance(ids, str):
            ids = [ids]
        userows = [2]
        # test if ids are genes to extract index
        if ids[0] in self.genes:
            for gene in ids:
                rows.append(self.genes.index(gene))
            userows.extend([x + self.skrows for x in rows])
            data = pd.read_csv(self.file, sep='\t', index_col=0, skiprows=lambda x: x not in userows)
            print(lambda x: x not in userows)
            data.drop(data.columns[range(0, self.hcols)], axis=1, inplace=True)
            return data, h
        # or if ids are wells to extract columns
        print(ids[0])
        if any([ids[0] in x for x in [h.index, h['addr'], h['well']]]):
            if ids[0] in h['well'].values:
                print('found match in well col ', ids[0])
                ids = [self.shortname + x for x in ids]
            if ids[0] in h['addr']:
                subh = h[h['addr'].isin(ids)]
            elif ids[0] in h.index.values:
                subh = h[h.index.isin(ids)]
            cols = subh.index.values
            for c in cols:
                coli.append(h.index.get_loc(c))
            coli = [x + (self.hcols+1) for x in coli]
            coli.insert(0, 0)
            data = pd.read_csv(self.file, sep='\t', index_col=0, skiprows=self.skrows+1,
                               usecols=coli, header=None)
            data.columns = cols
            print(data.columns)
            return data, subh


def main():
    path = '/Users/WRB/Desktop/PGA302-PGA303-PGA304_consensus.gct'
    d, h = dfsubset(path, ['200024_at'])
    d.head()

if __name__ == '__main__':
    main()

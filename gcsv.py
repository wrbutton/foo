#!/usr/bin/python
"""
12/13/16 - branched from gct into a csv class
                     will extract dataframe from the Median section of csv file
                     can also load a panel, and combine multiple plates for median/avg
                     saves a file by mimicing the structure of another csv file to
                     ensure pipeline compatibility
01/19/17 - added batch_adjust command to edit many csvs at once, and the
                     summarize_csvs to get basic metrics before processing

"""

import pandas as pd
import numpy as np
import collections as cll
import os, csv, gct, sys, pt, gt, math
import pickle


def cmbine_save(path, ctype):
    combined = get_combined(path, ctype)
    for filep in os.listdir(path):
        if filep.endswith('.csv'):
            c0 = Gcsv(filep)
            break
    csvf = os.path.join(path, c0.file)
    fname = '{}_{}.csv'.format(c0.shortname, ctype)
    output = os.path.join(path, fname)
    save_csv(csvf, combined, output)


def get_combined(path, ctype):
    # quick method to get median or mean of a panel
    panel = load_cpanel(path, pan=True)
    if ctype == 'median':
        combined = panel.median(axis=0)
        print('generating media of panel')
    elif ctype == 'mean':
        combined = panel.mean(axis=0)
        print('generating mean of panel')
    else:
        print('please enter valid combine type, returning blank')
        combined = []
    return combined


def load_cpanel(path, pan=False):
    # given folder path gets gct info and builds dframes for files
    # either merges into panel, or returns 2 dicts of gcts and dframes
    f_list = pt.get_flist(path, '.csv')
    gd, dfd = {}, {}
    for i, file in enumerate(f_list):
        gn, dfn = 'c' + str(i), 'df' + str(i)
        gd[gn] = Gcsv(file)
        dfd[dfn] = gd[gn].build_dframe()
    if pan == True:
        panel = pd.Panel({i: df for i, df in enumerate(dfd.values())})
        return panel
    elif pan == False:
        return gd, dfd


def summarize_csvs(path):
    """ provide path containing csv files to generate output summarizing levels 1 and 10
    for the plate as well as the posamp and ref """
    if path is None:
        path = gt.dflt_outpath(fldr_name='csv')
    results = cll.defaultdict(dict)
    f_list = pt.get_flist(path, '.csv')
    for file in f_list:
        try:
            c = Gcsv(file)
            d = c.build_dframe()
            results[c.shortname]['plate-L10'] = d['Analyte 10'].mean(axis=0)
            results[c.shortname]['Pos-L10'] = d.ix['B1']['Analyte 10']
            results[c.shortname]['Ref-L10'] = d.ix[['A2','B2']]['Analyte 10'].mean()
            results[c.shortname]['plate-L1'] = d['Analyte 1'].mean(axis=0)
        except:
            print('error with ' + file)
    res = pd.DataFrame(results)
    res = res.T
    outpath = os.path.join(path, 'csv_summary.txt')
    res.to_csv(outpath, sep='\t', float_format='%.0f')


# method to save a dafaframe object with gct file structure and headers
def save_csv(csvf, df, outpath):
    df = df.fillna('NaN')
    with open(csvf, 'rU') as in_file, open(outpath, 'w', newline='') as out_file:
        linereader = csv.reader(in_file, delimiter=',', quotechar='"')
        csvwriter = csv.writer(out_file, delimiter=',', quotechar='"',
                                                         quoting=csv.QUOTE_ALL)
        pd.options.display.float_format = '{i}'.format
        f = 0
        for line in linereader:
            # read through file until correct data section, writing uninteristing lines
            # to the corresponding edited file as-is
            if 'Data Type:' and 'Median' in line:
                csvwriter.writerow(line)
                # flag to read one more line and then change modes
                f = 1
            elif f == 1:
                csvwriter.writerow(line)
                break
            else:
                csvwriter.writerow(line)
        # read lines until blank checking wells and for matched ids from file
        for i in range(0,384):
            line = next(linereader)
            #dlist = ['{:.2f}'.format(x) if isnum(x) else 'NaN' for x in df.ix[i]]
            try:
                dlist = [x if gt.isnum(x) else 'NaN' for x in df.ix[i]]
                line[2:502] = dlist
                csvwriter.writerow(line)
            except IndexError:
                csvwriter.writerow(line)
                break
        # then just copy over the rest of the lines
        for line in linereader:
            csvwriter.writerow(line)


def check_11499(csvf):
    with open(csvf, 'rU') as in_file:
        linereader = csv.reader(in_file, delimiter=',', quotechar='"')
        pd.options.display.float_format = '{i}'.format
        for line in linereader:
            if line == ['DataType:', 'Count']:
                break
        line = next(linereader)
        line = next(linereader)
        c11 = line[12]
        c499 = line[500]
        fn = os.path.split(csvf)[-1]
        shn = '_'.join(fn.split('_')[0:2])
        print(fn, '\t', '11 count:', '{:4}'.format(c11), '\t', '499 count:', '{:4}'.format(c499))


def split_out_csv(csvf, d, outpath):
    """ takes in a larger csvfile to pull from, a dictionary 'd' mapping original to destination well,
    an outpath for the file. requires initial csv file, a dictionary of source: destination well mapping
    !!  that dictionary in 2 character well format !!
    counts are automatically scaled up 1.5x (but 499 and 11 left alone) """
    with open(csvf, 'rU') as in_file, open(outpath, 'w', newline='') as out_file:
        # set up csv file line writers and readers
        linereader = csv.reader(in_file, delimiter=',', quotechar='"')
        csvwriter = csv.writer(out_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_ALL)
        pd.options.display.float_format = '{i}'.format
        # initialize flags and variables
        q, f, csv_section = 0, 0, 0
        sample_names = []
        # define the csv file sections of interest
        sections = ['Median', 'Net MFI', 'Count', 'Avg Net MFI', 'Dilution Factor']
        # read through csv file until correct data section, writing uninteresting lines as is
        for line in linereader:
            # check for csv file sections of interest
            if any(sec in line for sec in sections) and 'DataType:' in line:
                # write that last line
                csvwriter.writerow(line)
                # flag to change file write mode to active data replacement
                f = 1
            elif f == 1:
                # write the header row just read
                csvwriter.writerow(line)
                # accept back the linereader object to continue moving through file
                # for line in linereader:
                for i in range(0, 384):
                    line = next(linereader)
                    # despite what pycharm says, needs to be line == [], NOT line is []
                    if line == []:
                        csvwriter.writerow(line)
                        break
                    # catch for the Avg Net MFI section, referencing the 'unknown2' names
                    if csv_section == 3:
                        if line[0] in sample_names:
                            csvwriter.writerow(line)
                            continue
                        else:
                            continue
                    # get the luminex well id
                    wid = line[0].split(',')[1].rstrip(')')
                    if wid in d.keys():
                        hdr = line[0].split(',')[0]
                        print(wid, '=>', d[wid])
                        nh = hdr + ',' + d[wid] + ')'
                        line[0] = nh
                        if csv_section == 2:
                            # store the 'unknown2' names to strip out unneeded for AVG MFI
                            sample_names.append(line[1])
                            # store 11 and 499 counts as is
                            c11, c499 = line[12], line[500]
                            # multiply counts to incorporate rptwls
                            line[2:502] = [int(float(x)) * 1.5 for x in line[2:502]]
                            # fill back in original counts
                            line[12], line[500] = c11, c499
                        csvwriter.writerow(line)
                # reset flag to be non-replace mode
                f = 0
                # increment the section header
                csv_section += 1
            else:
                csvwriter.writerow(line)
                continue


def adj_csv(csvf, outpath='dflt', mywells='all', pmap=None, adj_vect=None, separate=False, fixcounts=False,
            scalefact=1, adj_mode='center'):
    """ adjusts a csv file (whole or from list of wells) adjusting well values either by applying an adj_vector
    to each well, or by a scalefactor to each well. if a platemap is provided, then batch adjustment vectors will
    be applied. if separate is True only the wells in mywells list are written to file, otherwise only those are
    edited but all wells are copied over to new file. fixcounts ajdusts counts by 1.5, but if the file has been
    split out from a consolidated csv into sub RPTWLS files the count adjustment has been applied then

     new: adj_mode='center' is default and will apply adjustment vector (typically difference from medians),
     but if it's a list or tuples of numers you can combine a baseline subset value and plate value in a
     damping manner my specifying plate-based weight, local sample based-weight

     """
    if outpath is 'dflt':
        # outpath = csvf.strip('.csv') + '_adjusted.csv'
        if 'DP52' in csvf:
            if 'DP52_' in csvf:
                outpath = csvf.replace('DP52_', 'DP52_a')
            else:
                outpath = csvf.replace('DP52', 'DP52_a')
        elif 'DP53' in csvf:
            if 'DP53_' in csvf:
                outpath = csvf.replace('DP53_', 'DP53_a')
            else:
                outpath = csvf.replace('DP53', 'DP53_a')

    if mywells is 'all':
        mywells = get_2char_ids(pt.get_awells())
    else:
        mywells = get_2char_ids(mywells)
    if adj_vect is None:
        adj_vect = [0] * 500

    if len(mywells) > 34:
        print(print('adjusting ' + str(len(mywells))))
    else:
        print('adjusting ' + str(mywells))

    with open(csvf, 'rU') as in_file, open(outpath, 'w', newline='') as out_file:
        linereader = csv.reader(in_file, delimiter=',', quotechar='"')
        csvwriter = csv.writer(out_file, delimiter=',', quotechar='"',
                               quoting=csv.QUOTE_ALL)
        pd.options.display.float_format = '{i}'.format
        # initialize flags and variables
        flag, csv_section = 0, 0
        # define the csv file sections of interest
        hdr = ['Median', 'Net MFI', 'Count']
        # read through csv file until correct data section, writing uninteristing lines
        # to the corresponding edited file as-is
        for line in linereader:
            # check for csv file sections of interest
            if any(h in line for h in hdr):
                if 'DataType:' in line:
                    # write that last line
                    csvwriter.writerow(line)
                    # flag to change file write mode to active data replacement
                    flag = 1
            # pick up the next line and edit if it's a well of interest
            elif flag == 1:
                # write the header row just read, and loop through lines up to a full plate
                csvwriter.writerow(line)
                for i in range(0, 384):
                    # increment and get the next line
                    line = next(linereader)
                    # if a blank line (reached end of section) exit out
                    if line == []:
                        csvwriter.writerow(line)
                        break
                    # within MFI and Net MFI sections of csv
                    if csv_section < 2:
                        # want to only pull wells of interest to edit
                        well = line[0].split(',')[1].rstrip(')')
                        if well in mywells:
                            vals = [round(float(x) * scalefact) if x != 'NaN' else 0 for x in line[2:502]]
                            # if a plate map is provided, applies the per-batch adj-vector (which should be dict)
                            # otherwise applies a global adj_vector as a list to each well
                            if pmap is not None:
                                try:
                                    b = pmap[pmap['well'] == well]['batch'].values[0]
                                except IndexError:
                                    b = pmap[pmap['well'] == well]['batch'].values
                                z = zip(vals, list(adj_vect[b]))
                            else:
                                z = zip(vals, list(adj_vect))
                            # for regular adjustment vectors:
                            # combines scaled new with adj values using weights above
                            if adj_mode == 'center':
                                newvals = [abs(x[0] + x[1]) for x in z]
                            else:
                                try:
                                    p_wght = adj_mode[0]
                                    s_wght = adj_mode[1]
                                    print('adjusting plate to fill] in : ', p_wght, s_wght)
                                    newvals = [(p_wght * subv) + (s_wght * platev) for subv, platev in z]
                                except:
                                    print('if mode isnt center, it should be a weight between plate and subset')
                                    break
                                # insert new values into the line to be written
                            line[2:502] = newvals
                            csvwriter.writerow(line)
                            continue
                        elif well not in mywells:
                            if separate is False:
                                csvwriter.writerow(line)
                    # within the Count section of the csv, scale up if fixcounts is True
                    if csv_section == 2:
                        well = line[0].split(',')[1].rstrip(')')
                        if well in mywells:
                            if fixcounts is True:
                                # scale up counts, but store 11 and 499 counts as is
                                c11, c499 = line[12], line[500]
                                line[2:502] = [float(x) * 1.5 for x in line[2:502]]
                                line[12], line[500] = c11, c499
                                csvwriter.writerow(line)
                            elif fixcounts is False:
                                csvwriter.writerow(line)
                        elif well not in mywells:
                            if separate is False:
                                csvwriter.writerow(line)
                # reset flag to be non-replace mode
                flag = 0
                # increment the section header
                csv_section += 1
                continue
            else:
                csvwriter.writerow(line)


def get_2char_ids(wells):
    """ go from 3 char ids to the dumb luminex format """
    if isinstance(wells, str):
        wells = [wells]
    new_wells = []
    for w in wells:
        if len(w) == 3:
            if w[1] == '0':
                new_wells.append(w[0] + w[2])
            else:
                new_wells.append(w)
        elif len(w) == 2:
            new_wells.append(w)
    return new_wells


def get_3char_ids(wells):
    new_wells = []
    # wells = [w.split(',')[1].strip(')') for w in wells]
    for w in wells:
        if len(w) < 3:
            new_wells.append(w[0] + '0' + w[1])
        elif len(w) == 3:
            new_wells.append(w)
    return new_wells


def bulk_open_as_gct(path, drop_inv=False):
    # no current support for repeat wells files
    flist = pt.get_flist(path, '.csv')
    pdict = cll.defaultdict(list)
    for f in flist:
        shn = pt.get_shn(f)
        pdict[shn].append(f)
    for k, v in pdict.items():
        bsets = [x.split('_')[1] for x in v]
        if 'DP52' in bsets and 'DP53' in bsets and len(bsets) < 3:
            print(k, ' ok')
            df1 = open_as_gct(v[0], log=True)
            df2 = open_as_gct(v[1], log=True)
            df = pd.concat([df1, df2], axis=0)
            df.sort_index(inplace=True)
            if drop_inv is True:
                df = df[~df.index.str.contains('INV')]
            else:
                print('watch out, invariant genes included')
            outpath = os.path.join(path, k + '_pkexp.gct')
            gct.save_simplegct(df, outpath)
        else:
            print('error with: ', k)


def open_as_gct(file, log=False):
    csv = Gcsv(file)
    df = csv.transform_to_gct()
    #print(df.iloc[0,:])
    if log is True:
        # set 0's to 1s for log not to break
        try:
            df = df.drop('NaN')
        except:
            print('no NaN to drop')
        temp = df.copy()
        temp[temp==0] = np.nan
        minval = temp.min()[0]
        print(minval)
        df[df==0] = minval
        df = df.applymap(lambda x: round(math.log(x, 2), 3))
    return df


class Gcsv(object):
    # an object ossociated with extracting info from written gct files

    def __init__(self, infile):
        # pass the file argument as the filepath reference
        self.file = infile
        self.shortname = '_'.join(infile.split('/')[-1].split('_')[:2])

    def get_headers(self):
        # open file and extract the gct header information, storing dictionary values
        wells, hdrcl = [],[]
        with open(self.file, 'rU') as file:
            linereader = csv.reader(file, delimiter=',', quotechar='"')
            for i, line in enumerate(linereader):
                if 'Data Type:' and 'Median' in line:
                     self.skrows = i + 1
                     break
            for i, line in enumerate(linereader):
                if i >= 385:
                    self.samples = i
                    break
                if i == 0:
                    continue
                else:
                    if line == []:
                        self.samples = i
                        break
                    hdrcl.append(line[:2])
                    try:
                        well = line[0].split(',')[1].rstrip(')')
                        wells.append(well)
                    except IndexError:
                        print('')
                        print(line)
                        print('')
                        print('oops, no support yet for non-full csv files')
                        sys.exit()
            self.wells = wells
            self.hdrcl = hdrcl
            #print(self.wells)

    def null_empty(self, data):
        # used to replace process control wells with 'nan' so they won't plot
        for w in self.empty:
            data.iloc[:,w] = np.nan
        return data

    def build_dframe(self):
        # open file and extract the gct header information, storing dictionary values
        self.get_headers()
        # skip rows down to data, get first column value as index identifier
        data = pd.read_csv(self.file, sep=',', quotechar='"', skiprows=self.skrows,
                           nrows=self.samples-1, usecols=range(2,502))
        data.shortname = pt.get_shn(self.file)
        # retitle the index column name from Unnamed:0
        data = data.set_index([self.wells])
        # give back the completed dataframe with row + column ids
        return data

    def transform_to_gct(self):
        df = self.build_dframe()
        new_index = get_3char_ids(df.index)
        df.index = new_index
        # ref_file = '/Users/WRB/Dropbox/bin/python/analytes_to_probes_dict.p'
        ref_file = '/Users/WRB/Dropbox/bin/python/analytes_to_probes_dict.p'
        analyte_dict = pickle.load(open(ref_file, 'rb'))
        if 'DP52' in self.file:
            beadset_dict = analyte_dict['dp52']
        elif 'DP53' in self.file:
            beadset_dict = analyte_dict['dp53']
        genes = [col.lstrip('Analyte ') for col in df.columns]
        genes = [beadset_dict[int(analyte)] for analyte in genes]
        genes = ['0_' + x if 'INV' in x else x for x in genes]
        if 'DP52' in self.file:
            df['Analyte 499'] = 0
        elif 'DP53' in self.file:
            df['Analyte 11'] = 0
        df.columns = genes
        df.replace('NaN', 0)
        # del df['NaN']
        wells = get_3char_ids(df.index)
        df.index = wells
        df = df.T
        return df



#!/usr/bin/python

import pandas as pd
import numpy as np
import gct, sys, os, csv, math, glob, pickle, pt, gt, shutil, statistics
import collections as cll
import matplotlib.pyplot as plt 


def get_zscore(fpath, save=True):
    """ merged from separate zscore file. can either save the resulting file or return data
    the first fpath argument can be a file path or a [d, h] object already"""
    # basic setup
    if isinstance(fpath, str):
        g = gct.Gct(fpath)
        g.get_headers()
        df, h = gct.extractgct(fpath)
    else:
        try:
            df = fpath[0]
            h = fpath[1]
        except:
            print('error with path')

    zsd = cll.defaultdict(dict)

    for b in h['batch'].dropna().unique():
        if b == 'na':
          continue
        print('running zscore for {} batch {}'.format('pname', b))
        vw = gt.hsub(h, {'batch':b, 'type':'vehicle'}).index.values
        if len(vw) == 0:
            break
        veh = df[vw]
        # get median value across vehicle populations
        med = veh.median(axis=1)

        # populate the absolute deviation values per gene
        ad = cll.defaultdict(list)
        for v in veh.columns:
            for f in veh.index:
                ad[f].append(abs(med[f] - veh[v][f]))
        # assemble the median absolute value per gene
        mad = {}
        for k, v in ad.items():
            r = statistics.median(v)
            if 0 < r < 0.1:
                r = 0.1
            mad[k] = r
        # using the above progress though test and poscon wells
        # to calculate sample zscores
        tw = list(h[(h['batch'] == b) & (h['type'] == 'test')].index.values)
        pw = list(h[(h['batch'] == b) & (h['type'] == 'poscon')].index.values)
        wells = tw + pw
        for w in df[wells].columns:
            for feat in df.index:
                if mad[feat] == 0:
                    zs = 0
                else:
                    zs = (df[w][feat] - med[feat]) / (mad[feat] * 1.486)
                zsd[w][feat] = '{0:.3f}'.format(zs)

    # transform into dataframe, set index, null nonsense
    zsdf = pd.DataFrame(zsd)
    hs = h.loc[zsdf.columns]
    zsdf = zsdf.replace(['inf', '-inf'], np.nan).fillna('nan')
    if save is True:
        outpath = '{}_ZS.gct'.format(fpath.split('_', 1)[0])
        gct.save_headergct(zsdf, outpath, hs)
    else:
        return zsdf, hs


def get_dir_zscores(opath):
    # test whether input is file or directory, do single or loop through
    # all files as appropriate
    if os.path.isfile(opath) and opath.endswith('.gct'):
        get_zscore(opath)
    elif os.path.isdir(opath):
        for file in os.listdir(opath):
            if file[0] != '.':
                if file.endswith('.gct'):
                    targfile = os.path.join(opath, file)
                    get_zscore(targfile)


def create_rep_pairplots(df, h):
    """ assumes that the df + header contains all data of interest, preconcatenated multi plate ok """


def annotate_wells(path, wellfile, gene):
    """ reads from a text file ('MKD031:A03' one    per line) to grab in metadata
    from existing promiscuity and wellcorr files, and opens gct for zscore for given gene
    path for gct and annotation files"""
    # assemble file list, using list comp to remove rank gcts
    flg = pt.get_flist(path, '.gct')
    flg = [x for x in flg if 'rank' not in x]
    flg = [x for x in flg if 'consensus' not in x]
    # then build separate lists for the rest
    flr = pt.get_flist(path, 'ranks.gct')
    flc = pt.get_flist(path, 'consensus.gct')
    flp = pt.get_flist(path, 'promscty.txt')
    flwc = pt.get_flist(path, 'wellcorr.txt')
    data, h = pd.DataFrame(), pd.DataFrame()

    wd = cll.defaultdict(dict)

    hd = ['name', 'dose', 'rank', 'zs', 'wcc', 'pscr']

    with open(wellfile, 'rU') as f:
        for line in f:
            # selects the first well as key (name    dose    rank    *PGA103:A13*,PGA104:A15)
            wells = line.strip().split('\t')[-1].split(',')
            well = wells[0]
            plate = well.split(':')[0]
            print(well)
            # get first match of corresponding support files for the current plate
            gctfile = next(s for s in flg if plate in s)
            rankfile = next(s for s in flr if plate in s)
            confile = next(s for s in flc if plate in s)
            promfile = next(s for s in flp if plate in s)
            wellccfile = next(s for s in flwc if plate in s)
            # search through support files and extract info for given well
            # well replicate correlation
            for well in wells:
                with open(wellccfile, 'rU') as f1:
                    line = f1.readline().strip().split('\t')
                    #determine which column to search in
                    cccol = line.index('well cc')
                    name = line.index('name')
                    # detect if dose metadata present
                    try:
                        dose = line.index('dose')
                        doseflag = 1
                    except:
                        doseflag = 0
                    # loop through file, collecting well labels to see if on the list
                    for line in f1:
                        dat = line.strip().split('\t')
                        # if on the list, grab the correlation
                        if dat[0] == well:
                            # add well sub-dictionary parameter values
                            wd[well]['wcc'] = dat[cccol]
                            wd[well]['name'] = dat[name]
                            if doseflag == 1:
                                wd[well]['dose'] = dat[dose]
                            else:
                                wd[well]['dose'] = 'na'
                    try:
                        wd[well]['name']
                        break
                    except KeyError:
                        print('looking for alternate replicate in wellcc')
                        continue
            # next promiscuity
            with open(promfile, 'rU') as f2: 
                line = f2.readline().strip().split('\t')
                #determine which column to search in
                pcol = line.index('pscore')
                for line in f2:
                    dat = line.strip().split('\t')
                    if dat[0] == well:
                        wd[well]['pscr'] = dat[pcol]
            # on to zscore, grab from the consensus file
            g = gct.Gct(confile)
            g.get_headers()
            h = g.fh
            # use the except clause to use alternate replicate id in looking through
            # consensus files, in case it was used as pert header
            for w in wells:
                try:
                  wellcol = h.index.get_loc(w)
                  break
                except KeyError:
                  print('looking for alternate replicate in consensus')
                  continue
            coli = [wellcol + g.hcols + 1]
            coli.insert(0, 0)
            # grab just the index column and the sample column of interest as Series
            welldat = pd.read_csv(g.file, sep='\t', index_col=0, skiprows=g.skrows+1,
                                  usecols=coli, header=None)
            # grab the value in the Series for the given gene, add to dict
            wd[well]['zs'] = '{:.3f}'.format(welldat.ix[gene].values[0])
            # then do the same with ranks, being safe and re-interpreting new file info
            g = gct.Gct(rankfile)
            g.get_headers()
            h = g.fh
            coli = [wellcol + g.hcols + 1]
            coli.insert(0, 0)
            welldat = pd.read_csv(g.file, sep='\t', index_col=0, skiprows=g.skrows+1,
                                  usecols=coli, header=None)
            wd[well]['rank'] = welldat.ix[gene].values[0]

    # write out the assembled well dict values in order
    # name, dose, rank, zs, cc, prom, wells
    ofile = wellfile.replace('.txt', '-annot.txt')
    with open(wellfile, 'rU') as f:
        with open(ofile, 'w') as outf:
            wr = csv.writer(outf, delimiter='\t')
            for line in f:
                line = line.strip()
                wells = line.strip().split(',')
                for w in wells:
                    try:
                        d = wd[w]
                        row = []
                        for h in hd:
                            row.append(d[h])
                        row.append(line)
                        wr.writerow(row)
                        break
                    except KeyError:
                        continue


def get_reps(path, welllist, platelist, outpath):
    # platelist file required as reference input, same one as for consensus: d, PGA314, PGA315, PGA316
    # welllist is one well per line PGA103:14. code will look up in platelist and assemble
    # panel of those files header info to pull to search through for reps
    # output is rep wells list: PGA314:A14, PGA315:A15, PGA316:A18
    flist, fnames = pt.get_flist(path, 'ZS.gct', shn=True)
    prev_plates = ''
    with open(welllist, 'r') as infile, open(outpath, 'w') as outf:
        writer = csv.writer(outf)
        for line in infile:
            plates = []
            hdrs, reps = [], []
            wid = line.strip()
            plate = wid.split(':')[0]
            with open(platelist, 'r') as pl:
                for plateline in pl:
                    if plate in plateline:
                        plates = plateline.strip().split('\t')
                        break
                if plates == []:
                    print("didn't find plate in list: ", plate)
            #print(plates)
            # print(list(platelist))
            #plates = next(p for p in platelist if plate in p).strip().split('\t')
            print('plates: ', plates)
            #prev_plates = plates
            # if plates are same layout, finding reps is easy
            try:
                if len(plates) > 1:
                    if plates[0] == 's':
                        pids = plates[1:]
                        well = wid.split(':')[1]
                        reps = [p + ':' + well for p in pids]
                        writer.writerow(reps)
                    # if plates different layout, concat headers
                    elif plates[0] == 'd':
                        pids = plates[1:]
                        #if plates != prev_plates:
                        files = [fnames[p] for p in pids]
                        for f in files:
                            g = gct.Gct(f)
                            g.get_headers()
                            hdrs.append(g.fh)
                        h = pd.concat(hdrs, axis=0)
                        prev_plates = plates
                        #elif plates == prev_plates:
                        #    print('prev plates')
                        #    pass
                        # extract well ids based on match
                        # this could blow up if want to keep batches on each plate different
                        if 'dose' in h.columns:
                            n, d = h.loc[wid, ['name', 'dose']]
                            wset = h[(h['name'] == n) & (h['dose'] == d)].index.values
                        else:
                            n = h.loc[wid, 'name']
                            wset = h[h['name'] == n].index.values
                        writer.writerow(wset)
                    else:
                        print('> 1 but not d or s')
                # reverse of merge intraplate reps, calling get mergemode to determine
                # manner to combine replicates
                elif len(plates) == 1:
                    plate = plates[0].strip()
                    bmode, h = get_mergemode(fnames[plate])
                    if 'dose' in h.columns:
                        if bmode == 'acrossbatch':
                            mmode = 'plate_dose'
                            n, d = h.loc[wid, ['name', 'dose']]
                            wset = h[(h['name']==n) & (h['dose']==d)].index.values
                        elif bmode == 'withinbatch':
                            mmode = 'batch_dose'
                            n, d, b = h.loc[wid, ['name', 'dose', 'batch']]
                            wset = h[(h['batch'] == b) & (h['name'] == n) & (h['dose'] == d)].index.values
                        writer.writerow(wset)
                    elif 'dose' not in h.columns:
                        if bmode == 'acrossbatch':
                            mmode = 'plate'
                            n, d = h.loc[wid, ['name', 'dose']]
                            wset = h[h['name'] == n].index.values
                        elif bmode == 'withinbatch':
                            mmode = 'batch'
                            n, b = h.loc[wid, ['name', 'batch']]
                            wset = h[(h['batch'] == b) & (h['name'] == n)].index.values
                        writer.writerow(wset)
                    else:
                        print('< 1 but not dose')
            except KeyError:
                print(wid, ' well not found (failed well)')


def clean_zsfiles(path):
    for file in glob.iglob(path + '*ZSVCQNORM*'):
        nname = os.path.split(file)[1].split('_')[0] + '_ZS.gct'
        shutil.move(file, os.path.join(path, nname))
        print(nname)


def grab_rows(path, searchstring, ext, save=False):
    # extract rows from files in directory matching string
    flist = pt.get_flist(path, ext)
    outpath = os.path.join(path, searchstring + '_rows.txt')
    with open(outpath, 'w', newline='') as outf:
        for file in flist:
            with open(file, 'r') as f:
                for line in f:
                    if searchstring in line:
                        outf.write(os.path.split(file)[-1] + '\t')
                        outf.write(line)
                        

def well_consensus(wells, name):
    # takes list of pandas Series of Z-scored data as inputs
    # will return consensus of the minimum absolute score per gene
    consensus = []
    wl = len(wells)
    for g in wells[0].index:
        vals = [wells[i][g] for i in range(wl)]
        if all([v > 0 for v in vals]):
            consensus.append('{0:.3f}'.format(min(vals)))
        elif all([v < 0 for v in vals]):
            consensus.append('{0:.3f}'.format(max(vals)))
        else:
            consensus.append(0)
    consensus = pd.Series(consensus, index=wells[0].index, name=name)
    return consensus


def save_ccs(ccs, h, outpath):
    # handles the saving of replicate correlation information with header info
    # computes consensus 75th percentile of list of correlations to report
    # more compilcated as they're uneven lengths so can't dataframe it
    # !!! should round out the list of the ccs and do this by dataframe
    w = csv.writer(open(outpath, 'w'), delimiter='\t')
    headers = h.columns.values
    try:
        maxnumcc = max([len(x) for x in ccs.values()])
    except ValueError:
        print('no wells left, somethings amiss')
        return
    if maxnumcc > 5:
        maxnumcc = 5
    firstrow = ['well']
    firstrow.extend(h.columns.values)
    firstrow.append('well cc')
    if maxnumcc > 1:
        numresults = ['cc{}'.format(i+1) for i in range(maxnumcc)]
        firstrow.extend(numresults)
    w.writerow(firstrow)

    svnty = {}
    for k in ccs.keys():
        vnums = np.array([float(x) for x in ccs[k]])
        if len(vnums) == 0:
            svnty[k] = 0
        elif len(vnums) >= 1:
            svnty[k] = np.percentile(vnums, 75)
    skeys = sorted(ccs.keys(), key=lambda x: svnty[x], reverse=True)
    for key in skeys:
        row = [key]
        try:
            hdrs = h.loc[key].values
            row.extend(hdrs)
        except KeyError:
            #hdrs = h[h['well']==key][headers].values.ravel()
            hdrs = h[h['well']==key].ix[0].values 
            #hdrs = h[h['well']==key].values.ravel()
            row.extend(hdrs)
        row.append(svnty[key])
        if maxnumcc != 1:
            nums = [x for x in ccs[key]]
            if len(nums) > 5:
                nums = nums[:4]
            row.extend(nums)
        w.writerow(row)
    c = pd.Series(svnty)
    c.hist(range=(-.6,1))
    t = os.path.split(outpath)[-1].split('_')[0] + ' replicate cc'
    plt.title(t)
    plt.savefig(outpath.split('.')[0] + '.png')
    plt.clf()


def get_ccs(wells, ccs, wname, plot=True, outpath='dflt'):
    ccs[wname] = []
    for i in range(len(wells)):
        for j in range(1 + i, len(wells)):
            ccs[wname].append('{0:.2f}'.format(wells[i].corr(wells[j], method='pearson')))
    if plot is True:
        ax = sns.pairplot(wells)
        outpath = gt.dflt_outpath(outpath)
        plt.savefig(outpath + wname)
    return ccs


def wrapup_merge(new_wells, ccs, h, outpath):
    merged = pd.DataFrame(new_wells)
    gct.save_multiplateheadergct(merged, outpath, h)
    correlations = pd.Series(ccs).T
    correlations.sort_values(ascending=False, inplace=True)
    if '_consensus' in outpath:
        corr_out = outpath.replace('_consensus.gct', '_wellcorr.txt')
    else:
        corr_out = outpath.split('.')[0] + '_wellcorr.txt'
    save_ccs(ccs, h, corr_out)


def load_mergefiles(f):
    g = gct.Gct(f)
    df = g.build_dframe()
    h = g.fh
    new_wells = {}
    ccs = {}
    wids = []
    return g, df, h, new_wells, ccs, wids


def merge_wells(wids, df, new_wells, ccs):
    if len(wids) <= 1: 
        print('merge_wells 1 dropped failed well', wids)
        return new_wells, ccs
    wells, wid2 = [], []
    for w in wids:
        try:
            wells.append(df.loc[w])
        except KeyError:
            wid2.append(w)
            continue
    # should separate this out to a separate loop to try to grab the random
    # wells which somehow cannot be found, even though they exist
    wname = wids[0]
    if len(wells) <= 1: 
        print('merge_wells 2 dropped failed well', wells)
        return new_wells, ccs
    c = well_consensus(wells, wname)
    new_wells[wname] = c
    ccs = get_ccs(wells, ccs, wids[0])
    if len(wid2) > 0:
        new_wells, ccs = secondary_merge(wid2, df, new_wells, ccs)
    return new_wells, ccs


def secondary_merge(wids, df, new_wells, ccs):
    if len(wids) <= 1: 
        return new_wells, ccs
    wells = []
    for w in wids:
        try:
            wells.append(df.loc[w])
        except KeyError:
            print('shitting the bed on {}'.format(w))
            continue
    if len(wids) <= 1: 
        return new_wells, ccs
    wname = wids[0]
    c = well_consensus(wells, wname)
    new_wells[wname] = c
    ccs = get_ccs(wells, ccs, wids[0])
    return new_wells, ccs


def merge_intraplate_reps(f, shn, outpath, mmode):
    g, df, h, new_wells, ccs, wids = load_mergefiles(f)
    if mmode == 'plate_dose':
        for n in g.names:
            for d in h[h['name']==n]['dose'].dropna().unique():
                wset = h[(h['name']==n) & (h['dose']==d)]
                wids.append(wset.index.values)
    elif mmode == 'plate':
        for n in g.names:
            wset = h[h['name'] == n]
            wids.append(wset.index.values)
    elif mmode == 'batch':
        for b in g.batches:    
            for n in g.names:
                wset = h[(h['batch'] == b) & (h['name'] == n)]
                wids.append(wset.index.values)
    elif mmode == 'batch_dose':
        for b in g.batches:    
            for n in g.names:
                for d in h[h['name']==n]['dose'].unique():
                    wset = h[(h['batch'] == b) & (h['name'] == n) & (h['dose'] == d)]
                    wids.append(wset.index.values)
    df = df.T
    for wset in wids:
        wset = [wset[i] for i in range(len(wset))]
        new_wells, ccs = merge_wells(wset, df, new_wells, ccs)
    wrapup_merge(new_wells, ccs, h, outpath)


def swap_index(df, old, new):
    df[old] = df.index.values
    df.set_index(new, inplace=True)
    df[new] = df.index.values
    return df


def fill_out_headers(h, shn):
    if 'well' not in h.columns:
        if ':' in h.index.values[0]:
            h['well'] = h.index.str.split(':').str[1]
        else:
            h['well'] = h.index
    if 'addr' not in h.columns:
        if ':' in h.index.values[0]:
            h['addr'] = h.index.values    
        else:
            h['addr'] = shn + ':' + h.index.values
    return h


def merge_sameplate_reps(flist, shn, outpath):
    print('merging plate reps by well {}'.format(shn))
    # assume same well positions on all plates
    # load panel of gcts with common axes
    gd, dfd = gct.load_panel(flist, indv=True)
    headerlist = []
    platelist = []
    stuff = zip(gd.values(), dfd.values())
    wids, ccs, new_wells = [], {}, {}
    for g, df in stuff:
        h = g.fh
        h = fill_out_headers(h, g.shortname)
        headerlist.append(h)
        #if ':' not in df.columns[0]:
        #    df = pt.addr_id_df(df)
        platelist.append(df)
    df = pd.concat(platelist, axis=1)
    h = pd.concat(headerlist, axis=0)
    if ':' not in h.index.values[0]:
        print('switching')
        h = swap_index(h, 'well', 'addr')
    # get unique list of well ids contained in plate samples
    pwells = {wid.split(':')[1].strip() for wid in df.columns}
    pwells = sorted(pwells)
    for w in pwells:
        # select matches to well id
        wids = h[h['well']==w].index.values
        try:
            plate = h.loc[wids[0]]['addr'].split(':')[0]
        except IndexError:
            sys.exit()
        # plate = h[h['well']==w]['addr'].values[0].split(':')[0]
        subset = df[wids]
        wells = [subset[i] for i in subset]
        if len(wells) < 2:
            continue
        c = well_consensus(wells, w)
        wname = plate + ':' + w
        new_wells[wname] = c
        ccs = get_ccs(wells, ccs, wname)
    wrapup_merge(new_wells, ccs, h, outpath)


def merge_multiplate_reps(flist, outpath):
    print('merging multi plate reps by desc')
    # assume same well positions on all plates
    gd, dfd = gct.load_panel(flist, indv=True)
    headerlist = []
    platelist = []
    stuff = zip(gd.values(), dfd.values())
    for g, df in stuff:
        h = g.fh
        h = fill_out_headers(h, g.shortname)
        headerlist.append(h)
        # h = swap_index(h, 'well', 'addr')
        if ':' not in df.columns[0]:
            df = pt.addr_id_df(df)
        platelist.append(df)
    df = pd.concat(platelist, axis=1)
    h = pd.concat(headerlist, axis=0)
    wids, ccs, new_wells = [], {}, {}
    if 'dose' in h.columns:
        for n in h.name.unique():
            for d in h[h['name']==n]['dose'].unique():
                wset = h[(h['name'] == n) & (h['dose'] == d)]['addr'].values
                wids.append(wset)
    else:
        for n in h.name.unique():
            wset = h[h.name == n]['addr'].values
            wids.append(wset)
    df = df.T
    df = df.sort_index()
    for wset in wids:
        new_wells, ccs = merge_wells(wset, df, new_wells, ccs)
    wrapup_merge(new_wells, ccs, h, outpath)


def get_mergemode(f):
    g = gct.Gct(f)
    g.get_headers()
    h = g.fh
    numreps = []
    if 'dose' in h.columns:
        uniques = []
        for n in g.names:
            for d in h[h['name']==n]['dose'].unique():
                uniques.append((n, d))
        for b in g.batches:
            for u in uniques:
                n, d = u
                num = len(h[(h['batch']==b) & (h['name']==n) & (h['dose']==d)].index)
                numreps.append(num)
        if sum(numreps)/len(numreps) > 1.5:
            bmode = 'withinbatch'
        else:
            bmode = 'acrossbatch'
    else:
        for b in g.batches:
            for n in g.names:
                num = len(h[(h['batch']==b) & (h['name']==n)].index)
                numreps.append(num)
        if sum(numreps)/len(numreps) > 1.5:
            bmode = 'withinbatch'
        else:
            bmode = 'acrossbatch'
    return bmode, h


def bulk_merge(path, mergelist):
    flist, shnd = pt.get_flist(path, '.gct', shn=True)
    platereps = pt.txt2list(mergelist)
    for pset in platereps:
        print(pset)
        if len(pset) > 1:
            outpath = os.path.join(path, '-'.join(pset[1:]) + '_consensus.gct')
            files = [shnd[n] for n in pset[1:]]
            if pset[0] == 's':
                pset = pset[1:]
                merge_sameplate_reps(files, pset,    outpath)
            elif pset[0] == 'd':
                pset = pset[1:]
                merge_multiplate_reps(files, outpath)
        else:
            try:
                f = shnd[pset[0]]
            except KeyError:
                print('file {} not found, skipping'.format(pset[0]))
                continue
            shn = [pset[0]]
            bmode, h = get_mergemode(f)
            outpath = os.path.join(path, pset[0] + '_consensus.gct')
            if 'dose' in h.columns:
                if bmode == 'acrossbatch':
                    mmode = 'plate_dose'
                elif bmode == 'withinbatch':
                    mmode = 'batch_dose'
            else: 
                if bmode == 'acrossbatch':
                    mmode = 'plate'
                elif bmode == 'withinbatch':
                    mmode = 'batch'
            print(mmode)
            merge_intraplate_reps(f, shn, outpath, mmode)


def get_promiscuity(df, t=0):
    """ sum absolute val of zscores above threshold in dataframe, returns dictionary """
    pscores = cll.OrderedDict()
    df = abs(df)
    for w in df.columns:
        pscores[w] = int(round(sum([x for x in df[w] if x >= t])))
    return pscores


def survey_promiscuity(path, t=0, indv=True):
    """ goes through a directory, if indv=True saves per-well promiscuity for each file, but
    also then saves a composite summary of range of percentiles across plates """
    flist = pt.get_flist(path, '.gct')
    prom_summary = {}
    quantiles = [1, .95, .9, .75, .5, .25, .1, .05, .01]
    for file in flist:
        d, h = gct.extractgct(file)
        prom_scores = get_promiscuity(d, t)
        f = pt.get_shn(file)
        if indv is True:
            ps = pd.DataFrame(prom_scores, index=[0]).T
            ps.columns = ['promiscuity']
            ps = pd.merge(h, ps, left_on='addr', right_index=True, how='outer')
            outpath = os.path.join(path, f + '_promiscuity.txt')
            ps.to_csv(outpath, sep='\t')
        # otherwise continue with overall summary
        l = np.array([int(x) for x in prom_scores.values()])
        prom_summary[f] = [np.percentile(l, x*100) for x in quantiles]
    pdf = pd.DataFrame(prom_summary, index=quantiles)
    outpath = os.path.join(path.split('.')[0], '_promiscuity_summary.txt')
    pdf.transpose().to_csv(outpath, sep='\t')


def get_gene_activity(df, threshs):
    activity = cll.OrderedDict()
    for t in threshs:
        val = df[abs(df) > t].count().sum()
        activity[t] = val
    return activity


def gene_activities(df, outpath):
    t = [2, 4, 7, 10, 15]
    adf = pd.DataFrame(columns=t)
    m = np.array([1, 5, 10, 15, 20])
    for g in df.index:
        activity = get_gene_activity(df.ix[g], t)
        adf.loc[g] = activity
    adf['actv score'] = (adf * m).sum(axis=1)
    adf = adf.sort_values(['actv score'], ascending=False)
    adf[adf < 1] = np.nan
    adf.replace(np.nan, '', inplace=True)
    adf.to_csv(outpath, sep='\t')


def survey_activities(path):
    # set the promiscuity score threshold to filter out in quiet activity
    pthresh = 1000
    # flist = glob.iglob(path + '*consensus.gct')
    flist = pt.get_flist(path, '.gct')
    for f in flist:
        print('getting gene activity for {}'.format(os.path.split(f)[-1]))
        g = gct.Gct(f)
        d = g.build_dframe()
        outpath = os.path.join(path, g.shortname + '_gene_actvty.txt')
        gene_activities(d, outpath)
        outpath = os.path.join(path, g.shortname + '_gene_qactvty.txt')
        pscores = get_promiscuity(d)
        pdf = pd.DataFrame(list(pscores.items()), columns=['well', 'pscore'])
        pdf.set_index(['well'], inplace=True)
        if 'dose' in g.fh.columns:
            h = g.fh[['name', 'dose', 'batch']]
        else:
            h = g.fh[['name', 'batch']]
        pfh = pd.concat([h, pdf], axis=1, join='inner')
        pfh.sort_values(['pscore'], ascending=False, inplace=True)
        pout = os.path.join(path, g.shortname + '_promscty.txt')
        pfh.to_csv(pout, sep='\t')
        qwells = [w for w in pscores if pscores[w] < pthresh]
        qdf = d[qwells]
        gene_activities(qdf, outpath)


def summarize_act(file):
    r = pd.read_csv(file, sep='\t', header=None)
    # create header for file
    rc = ['plate','gene']+['zs-'+str(x) for x in [2,4,7,10,15]]+['actv score']
    r.columns = rc


def rank_files(path):
    """ reads a zscore gct file (probably consensus) and saves file of ranks
            of each gene within each sample - from highest to lowest values"""
    flist = pt.get_flist(path, '.gct')
    for f in flist:
        g = gct.Gct(f)
        print(f)
        d = g.build_dframe()
        ds = d.rank(ascending=False)
        outpath = os.path.join(os.path.split(path)[0], g.shortname + '-ranks.gct')
        gct.save_headergct(ds, outpath, g.fh)


def grab_ranks(path, feat, hilo=1, t=40):
    """ survey folder for the wells in whic given gene istop ranked wells of given feat by sorted z score). generates 
        overall list ofdefault rank output is in descending order, highest zscore = 1
        hilo: 1 = high, upregulated genes (default rank order)
                    0 = low, downnregulated genes """
    outpath = os.path.join(path, '_rank_summary.txt')
    flist = pt.get_flist(path, 'ranks.gct')
    # set dummy starting point for low rank
    lowest = 500
    # create blank template dataframe
    summary = pd.DataFrame()        
    for f in flist:
        d, h = gct.extractgct(f)
        # flip rank order as needed
        if hilo > 1:
            d = 978 - d
        # get column ids for ranks below threshold
        wells = d.columns[d.ix[feat] < t]
        # extract portion of dataframe
        ranks = d.ix[feat, wells]
        ranks = pd.DataFrame(ranks)
        # assign plate column to each well id entry, re-order cols
        ranks['plate'] = pt.get_shn(f).split('-')[0]
        # concat portion to overall dataframe
        summary = pd.concat([summary, ranks])
        # check and store the lowest rank
        newlow = min(d.ix[feat])
        if newlow < lowest:
            lowest = newlow
    # re-shuffle the column order
    summary['well'] = summary.index.values
    summary = summary[['plate', 'well', feat]]
    print('\n', feat, int(lowest))
    summary.to_csv(outpath, sep='\t', index=None)


def get_gene_modulators(path, plist, gene, gthresh, pthresh=None):
    # searches for wells which modulate the given gene above absolute zs threshold
    # extracts well correlation and promiscuity files in addition to zs of gene
    # 1. just get a list of the plates i care about in the directory
    shnamelist = pt.txt2list(plist)
    # then strip them, by taking first element
    shnamelist = [x[0] for x in shnamelist]
    files = glob.glob(path + '*.gct')
    outpath = os.path.join(path, '_gene_modulators.txt')
    with open(outpath, 'w') as outf:
        wr = csv.writer(outf, delimiter='\t')
        wr.writerow(['plate', 'well', 'name', 'dose', 'batch', 'zs', 'well cc', 'promiscuity'])
        for f in files:
            fshn = pt.get_shn(f)
            if fshn in shnamelist:
                print('ok', f)
                # grab the wells which modulate gene above threshold, and their zs val
                d, h = gct.extractgct(f)
                wells = d.columns[abs(d.loc[gene]) > gthresh].values
                print(wells)
            else:
                continue
            # assemble wells and vals into dictionary
            wd = cll.defaultdict(dict)
            # populate the dictionary with info: winfo['A03']['zs'] = -4.3 (zsval in A03)
            for w in wells:
                wd[w]['zs'] = d.ix[gene, w]
            # then go collect corr and promiscuity info
            f1n = fshn + '_wellcorr.txt'
            f2n = fshn + '_promscty.txt'
            fpath1 = os.path.join(path, f1n)
            with open(fpath1, 'rU') as f1: 
                line = f1.readline().strip().split('\t')
                #determine which column to search in
                col = line.index('well cc')
                # !!! is this only for wells vs address???
                # loop through file, collecting well labels to see if on the list
                for line in f1:
                    dat = line.strip().split('\t')
                    # if on the list, grab the correlation
                    if dat[0] in wells:
                        wd[dat[0]]['wcc'] = dat[col]
            fpath2 = os.path.join(path, f2n)
            with open(fpath2, 'rU') as f2: 
                line = f2.readline().strip().split('\t')
                #determine which column to search in
                col = line.index('pscore')
                dd = {}
                desc = ['name', 'dose', 'batch']
                for n in desc:
                    if n in line:
                        dd[n] = line.index(n)
                for line in f2:
                    dat = line.strip().split('\t')
                    if dat[0] in wells:
                        wd[dat[0]]['pscr'] = dat[col]
                        for d in desc:
                            try:
                                wd[dat[0]][d] = dat[dd[d]]
                            except:
                                wd[dat[0]][d] = 'na'
            # write one row per plate:well match
            for w in wd.keys():
                row = [fshn, w]
                hdr = ['name', 'dose', 'batch', 'zs', 'wcc', 'pscr']
                # if pthresh defined, skip those vals
                if pthresh != None:
                    if int(wd[w]['pscr']) > pthresh:
                        print('too loud, skipping')
                        continue
                for n in hdr:
                    row.append(wd[w][n])
                wr.writerow(row)




def main():

    mergelist = '/Users/wrb/Desktop/GJB_merge/mergelist.txt'
    #mergelist = '/Volumes/WRBHDD/wrb/Desktop/foolist.txt'

    #path = '/Volumes/WRBHDD/wrb/Desktop/processed/'
    path = '/Users/wrb/Desktop/GJB_merge/'
 
    plist = '/Volumes/QCexternal/WRBfiles/haplo/plates.txt'
    #plist = '/Volumes/WRBHDD/wrb/Desktop/platelist.txt'

    ''' first run bulk merge on ZS files in folder, along with merge key file
            listing either one plate per line or 's' or 'd' along with multiple 
            plates (tab sep between all), for same trt layout or different '''
    bulk_merge(path, mergelist)
    
    ''' then remove individual ZS gcts from folder, leaving only consensus
            survey activities across the consunsus gct files, which produces the
            activity and qactivity summaries, along with promiscuity '''
    #survey_activities(path)

    ''' grab ranks se'''
    # rank_files(path)

    ''' then grab rows to assess the behavior of the target gene of interest
            across plates, looking at the activity vs qactivity measures.
            decide on which plates are of subsequent interest '''
    #grab_rows(path, '200678_x_at')

    ''' once a range of plates has been narrowed down from which you want to
            pull individual wells from, then create text file listing plates of interest.
            Then run code with gene and (absolute) zscore threshold to extract mods of.
            The output will grab the promiscuity score, well cc, zscore of each well. '''
    #get_gene_modulators(path, plist, '200678_x_at', 3)

    '''
    then figure out what you care about in terms of wells from there, need to then 
    run them through the 'gather reps' thing to assemble the sets of individual

    '''
if __name__ == '__main__':
    main()
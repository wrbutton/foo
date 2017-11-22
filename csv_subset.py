
import gt, gcsv, os, csv, sys
import pandas as pd


def split_rename_csvs(path='dflt', mapfile='dflt'):
    """ updated October 2017, relies upon other gcsv functions
    but all in one batch and split them apart into individual csv files which can
    be run as RPTWLS. Requires input excel file of form below, can list many files.
    
    assumes to be in target directory named 'mapping.xlsx'
    Can use same Det Well as Destination well if simply want to split out
    
    CSV Name	            Det Well	Destination
    ITP549_DP52_RPTWLS.csv	E1	                C11
    ITP549_DP52_RPTWLS.csv	F1	                D11
    ITP553_DP52_RPTWLS.csv	I1	                C11
    ITP553_DP52_RPTWLS.csv	J1	                D11
    """
    if path is 'dflt':
        path = '/Users/WRB/Desktop/splitcsv/'
    if mapfile is 'dflt':
        mapfile = gt.get_flist(path, 'mapping.xlsx')[0]

    try:
        # load the plate/well remapping file
        m = pd.read_excel(mapfile)
        dest_csvs = m['CSV Name'].unique()
        src_csvs = gt.get_flist(path, '.csv')
    except:
        sys.exit(" dflt path is '/Users/WRB/Desktop/splitcsv/' and assumes 1 csv and mapping.xlsx")

    try:
        src_52_csv = [x for x in src_csvs if 'DP52' in x and 'RPTWLS.csv' not in x][0]
    except:
        print('no d52 source csv')
    try:
        src_53_csv = [x for x in src_csvs if 'DP53' in x and 'RPTWLS.csv' not in x][0]
    except:
        print('no d53 source csv')

    # loop through destination csvs pulled from mapping file
    for csv_file in dest_csvs:
        # get subset of the plate remapping dataframe for the same plate
        subm = m[m['CSV Name']==csv_file]
        outpath = os.path.join(path, csv_file)
        # convert dataframe of excel remapping into an easy dictionary
        dest = list(subm['Destination'].values)
        src = gcsv.get_2char_ids(list(subm['Det Well'].values))
        dest = gcsv.get_2char_ids(list(subm['Destination'].values))
        mydict = dict(zip(src, dest))
        # then call on the split_csv command in gcsv module
        if 'DP52' in csv_file:
            print(os.path.split(src_52_csv)[-1] + ' => ' + csv_file)
            gcsv.split_out_csv(src_52_csv, mydict, outpath)
        elif 'DP53' in csv_file:
            print(os.path.split(src_53_csv)[-1] + ' => ' + csv_file)
            gcsv.split_out_csv(src_53_csv, mydict, outpath)
        else:
            print('no match')


def adj_rptwls_to_plates(platepath='dflt', subpath='dflt'):
    """ batch adjust individual repeat wells files, assuming single batch per plate
    derives adjustment vector for each RPTWLS compared to full plate for bead set and runs
    if there are plates that have multiple batches, the plate map for those should be placed
    in the plates folder, and if present adj_vectors will be made on a per-batch basis """

    if platepath is 'dflt':
        platepath = '/Users/WRB/Desktop/plates/'
    if subpath is 'dflt':
        subpath = '/Users/WRB/Desktop/rpts/'

    pflist = gt.get_flist(platepath, '.csv')
    sflist = gt.get_flist(subpath, '.csv')
    sflist = [x for x in sflist if 'adjusted' not in x]
    maplist = [x for x in gt.get_flist(platepath, '.xlsx')]

    pshn, sshn = {}, {}

    for p in pflist:
        sh = '_'.join(os.path.split(p)[1].split('_')[:2])
        pshn[sh] = p
    for s in sflist:
        sh = '_'.join(os.path.split(s)[1].split('_')[:2])
        sshn[sh] = s

    for s in set(sshn.keys()):
        print(s)
        platefile = pshn[s]
        subfile = sshn[s]
        try:
            mapfile = [x for x in maplist if s.split('_')[0] in x][0]
        except:
            mapfile = None

        if mapfile:
            print('adjusting ' + s + ' by batch with map')
            m = pd.read_excel(mapfile)
            batches = m['batch'].dropna().unique()
            badjvect = {}
            dp = gcsv.open_as_gct(platefile)
            ds = gcsv.open_as_gct(subfile)
            subwells = ds.columns.values
            pwells = [x for x in dp.columns.values if x not in subwells]
            for b in batches:
                wids = m[m['batch'] == b]['well'].values
                pids = [x for x in wids if x in pwells]
                sids = [x for x in wids if x in subwells]
                pbmed = dp[pids].median(axis=1)
                sbmed = ds[sids].median(axis=1)
                badjvect[b] = round(pbmed - sbmed)
            gcsv.adj_csv(subfile, adj_vect=badjvect, pmap=m)
        # if no mapfile present, proceed with single batch adj_vect
        else:
            print('adjusting ' + s)
            d = gcsv.open_as_gct(platefile)
            pmed = d.median(axis=1)
            d = gcsv.open_as_gct(subfile)
            submed = d.median(axis=1)
            diff_vect = round(pmed - submed)
            gcsv.adj_csv(subfile, adj_vect=diff_vect)


def adj_fillin_wells(path='dflt', wells='top,bot'):
    """ adjust full csv files for the traditional Turbocapture problem areas, default top and bottom
     if a given plate map file is present in the folder for a given plate, then individual batch
     adj_vectors are prepared and passed through"""

    if path is 'dflt':
        path = '/Users/WRB/Desktop/plates/'

    try:
        pflist = gt.get_flist(path, '.csv')
        maplist = gt.get_flist(path, '.xlsx')
        pflist = [x for x in pflist if 'adjusted' not in x]
    except:
        sys.exit("dflt path is '/Users/WRB/Desktop/plates/' and top/bot = True")

    mywells = []
    if 'top' in wells:
        mywells.extend(gt.well_range("C11", "F14"))
    if 'bot' in wells:
        mywells.extend(gt.well_range("K11", "N14"))
    else:
        mywells = wells

    mshn = {}
    for m in maplist:
        sh = gt.get_shn(m).strip('.xlsx')
        print(sh)
        mshn[sh] = m
    print(maplist)

    for platefile in pflist:
        plate_name = gt.get_shn(platefile)
        d = gcsv.open_as_gct(platefile)
        # if excel map file is in directory, cycle through to get batch adj vectors
        if plate_name in mshn.keys():
            print('adjusting ' + plate_name + ' by batch with map')
            m = pd.read_excel(mshn[plate_name])
            batches = m['batch'].dropna().unique()
            badjvect = {}
            pwells = [x for x in d.columns.values if x not in mywells]
            for b in batches:
                wids = m[m['batch'] == b]['well'].values
                pids = [x for x in wids if x in pwells]
                sids = [x for x in wids if x in mywells]
                pbmed = d[pids].median(axis=1)
                sbmed = d[sids].median(axis=1)
                badjvect[b] = round(pbmed - sbmed)
            gcsv.adj_csv(platefile, mywells=mywells, adj_vect=badjvect, pmap=m)
        # otherwise just proceed with single batch adjustment
        else:
            print('adjusting ' + plate_name)
            pmed = d.drop(mywells, axis=1).median(axis=1)
            sub = d[mywells]
            submed = sub.median(axis=1)
            diff_vect = round(pmed - submed)
            gcsv.adj_csv(platefile, mywells=mywells, adj_vect=diff_vect)




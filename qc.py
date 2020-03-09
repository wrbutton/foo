
from PIL import Image
import os, glob, gt, shutil, errno, boto3, subprocess, gct, csv
import pandas as pd
from sklearn import decomposition
import matplotlib.pyplot as plt
import collections as cll


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


def check_final(path='dflt'):
    """ check numbers of row/columns, number of fails and decimal places of final data files """
    if path == 'dflt':
        path = gt.dflt_outpath(fldr_name='finaldata')

    f_list = gt.get_flist(path, ext='.gct')

    for file in f_list:
        g = gct.Gct(file)
        g.get_headers()
        try:
            txt = g.file.split('.')[0] + '.txt'
        except:
            try:
                txt = g.file.split('.')[0] + '.xlsx'
            except:
                pass

        try:
            print(sub_check_failed(g, txt))
            fails, fail_result = sub_check_failed(g, txt)
            result = sub_check_lines(g) and sub_check_columns(g) and fail_result
            dplaces = sub_check_decimal(g)
        except FileNotFoundError:
            result = False
            fails = 'no map!!'

        print('{} - {} - {} failed wells - {} dplaces'.format(g.shortname, result, fails, dplaces))


def check_data(path='dflt'):
    """ a better final map checker """
    if path == 'dflt':
        path = gt.dflt_outpath(fldr_name='finaldata')

    flist = gt.get_flist(path, ext='.gct')
    maplist = gt.get_flist(path, ext='.txt')
    maplist.extend((gt.get_flist(path, ext='.xlsx')))

    for f in flist:
        shn = gt.get_shn(f).split('.')[0]
        try:
            mapfile = [x for x in maplist if shn in x][0]
        except:
            print(f'error with map file {shn}')

        g = gct.Gct(f)
        g.get_headers()
        g.get_wells()
        datwells = g.wells

        mymap = gct.extractmap(mapfile)

        mapwells = gt.hsub(mymap, {'type':['vehicle', 'poscon', 'test']})['well'].values

        res = set(datwells) - set(mapwells)

        if len(res) == 0:
            print(f'{shn} ok, {380-len(datwells)} failed wells')
        else:
            print(f'eror with map/data {shn}, {len(datwells)}/{len(mapwells)}')


def combine_fails(path='dflt', ret=False, summ=False, sep=False, thresh=1):
    if path == 'dflt':
        path = gt.dflt_outpath(fldr_name='newQC')
    fl = gt.globit(path, '*QC_fail*')
    files = ' '.join(fl)
    #cmd_str = 'cat ' + files + ' > ' + os.path.join(path, 'QC_fail.txt')
    #subprocess.run(cmd_str, shell=True)
    datlist = []
    for f in fl:
        dat = pd.read_csv(f, sep='\t', skiprows=1)
        dropcols = [x for x in dat.columns if 'Unnamed' in x]
        dat = dat.drop(dropcols, axis=1)
        dat.dropna(inplace=True)
        try:
            dat = dat[dat['Batch'] != ' ']
        except:
            pass
        if sep == False:
            try:
                dat = dat[dat['Batch'] != 'Batch']
            except:
                pass
        datlist. append(dat)
    data = pd.concat(datlist, axis=0)
    data.to_csv(os.path.join(path, 'QCfail_summary.txt'), sep='\t')

    if summ is True:
        gbname = data.groupby('PERT_DESC').size()
        print(gbname[gbname > thresh])
        gbbatch = data.groupby('Batch').size()
        print(gbbatch[gbbatch > thresh])

        # this subsets down to show how many doses totally fail (3 reps each) per name
        # g = f.groupby(['PERT_DESC', 'DOSE']).size()
        # res = g[g > 2].groupby('PERT_DESC').size().sort_values(ascending=False)

    if ret is True:
        return data


def assemble_ref_dat(path):
    """ to gather together all reference RNA wells within the given path """
    fl = gt.globit(path, '*_ref_n*')
    dl, hl = [], []
    for f in fl:
        dr, hr = gct.extractgct(f)
        dr, hr = gt.dsub(dr, hr, {'well':['A02','B02']})
        dr = round(dr, 2)
        dl.append(dr)
        hl.append(hr)
    alldata = pd.concat(dl, axis=1)
    return alldata


def run_ref_pca(path='dflt', data='dflt', save=True, label=True):
    """ requires pointing to a pre-assembled data file to reference new pca against"""
    if path == 'dflt':
        path = gt.dflt_outpath(fldr_name='newQC')
    refdf = assemble_ref_dat(path)
    ns = len(refdf.columns)
    if isinstance(data, str):
        if data == 'dflt':
            dat_file = '/Users/WRB/Dropbox/Areas of Focus/_Genometry/Analysis Projects/reference_ref_pca.csv'
            bulk_data = pd.read_csv(dat_file, delimiter='\t', index_col=0)
    else:
        bulk_data = data
    try:
        newcols = [x for x in refdf.columns if 'B01' not in x]
        refdf = refdf[newcols]
    except:
        pass
    alldata = pd.concat([bulk_data, refdf], axis=1)
    pca = decomposition.PCA(n_components=2).fit_transform(alldata.T)
    fig, ax = plt.subplots()
    ax.scatter(pca[:, 0], pca[:, 1], c='grey')
    ax.scatter(pca[-ns:, 0], pca[-ns:, 1], c='red')

    if label is True:
        for l, x, y in zip(refdf.columns.values, pca[-ns:, 0], pca[-ns:, 1]):
            plt.annotate(l, xy=(x, y), xytext=(-2,2),
                         textcoords='offset points', ha='right', va='bottom')

    ax.set_title('REF RNA PCA')
    #pca = pd.DataFrame(pca, index=alldata.columns, columns=['comp1','comp2'])
    if save is True:
        outpath = gt.dflt_outpath(fldr_name=None, fn='pca.png')
        plt.savefig(outpath)


def check_s3(p):
    """ use aws cli in unix to check contents of s3 folder, p = path. p = 'd' is dflt '_review' """
    if p == 'd':
        p = '_review'
    cmd_str = f'aws s3 ls s3://genometry/{p}/'
    subprocess.run(cmd_str, shell=True)


def dl_data(sc='q', src='dflt', dest='dflt', search=None, excl=None, ext=None, p=False):
    """ download data from s3. 'sc' is shortcut, can be 'q' for qc, 'g' for gct, 'z' for zscore,
    'e' for enrichment, 'f' for final """
    if dest == 'dflt':
        dest = gt.dflt_outpath(fldr_name='newQC')
    elif dest == 'foo':
        dest = gt.dflt_outpath(fldr_name='foo')
    else:
        if '/' in dest or '\\' in dest:
            try:
                os.mkdir(dest)
            except:
                pass
        else:
             dest = gt.dflt_outpath(fldr_name=dest)

    tempdest = gt.dflt_outpath(fldr_name='temp_copy_transfer')

    s3c = boto3.client('s3')
    pref = '_review/'
    if src == 'dflt':
        items = s3c.list_objects(Bucket='genometry', Prefix=pref, Delimiter='/')
        folds = sorted([list(x.values())[0].replace(pref, '').strip('/') for x in items['CommonPrefixes']])
        if len(folds) == 0:
            print('hm, zero files in list')
        if len(folds) == 1:
            fold = folds[0]
        if len(folds) > 1:
            fold = folds[-1]
        print('downloading from ', fold)
        src = 's3://genometry/' + pref + fold
        # grab PCA ppt and latest coordinates txt
        s3c.download_file('genometry', 'PCA_analysis/PCA2.pptx', os.path.join(dest, 'PCA.pptx'))
        coords = s3c.list_objects(Bucket='genometry', Prefix='PCA_analysis/')
        coord = sorted([x['Key'] for x in coords['Contents']])[-1]
        s3c.download_file('genometry', coord, os.path.join(tempdest, 'PCA_coords.txt'))
    else:
        src = 's3://genometry/' + src

    search_args = []

    # parse shortcut
    if 'q' in sc:
        search_args.append(('*_qc/*', ''))
    if 'g' in sc:
        search_args.append((['*_fullqnorm_*', '*_QNORM_sorted*', '*_ref_n*', '*_Qctrl_n*'], ''))
        ext = '.gct'
    if 'z' in sc:
        search_args.append(('*_ZSVCQNORM_*', ''))
        ext = '.gct'
    if 'e' in sc:
        search_args.append(('*_escore/*', '*.gct'))
    if 'f' in sc:
        search_args.append(('*_final/*', ''))
        dest = gt.dflt_outpath(fldr_name='finaldata')
        ext = ['.gct', '.txt']
    if 'i' in sc:
        search_args.append(('*_ZSVCINF_*', ''))
        ext = '.gct'
    search_args.append(('', excl))

    if search is not None:
        if '*' not in search:
            search = '*' + search + '*'
        search_args.append((search, ''))

    if excl is not None:
        if '*' not in excl:
            excl = '*' + excl + '*'
        search_args.append(('', excl))

    for search, excl in search_args:
        cmd_str = f'aws s3 cp --recursive {src} {tempdest} --exclude "*"'
        if isinstance(search, str):
            search = [search]
        for st in search:
            cmd_str += f' --include {st}'
        if excl != '':
            cmd_str += f' --exclude {excl}'

        print(cmd_str)

    subprocess.run(cmd_str, shell=True)

    if ext is not None:
        if isinstance(ext, str):
            ext = [ext]
        for ex in ext:
            fl = gt.globit(tempdest, f'*{ex}')
            for f in fl:
                file_dest = os.path.join(tempdest, os.path.basename(f))
                shutil.move(f, file_dest)
            subdirs = [x[0] for x in os.walk(tempdest)][1:]
        for sd in subdirs:
            try:
                shutil.rmtree(sd)
            except:
                pass
    # do final copy from temp to main destination, erase temp
    for f in gt.get_flist(tempdest):
        file_dest = os.path.join(dest, os.path.basename(f))
        shutil.move(f, file_dest)
    shutil.rmtree(tempdest)

    if p is True:
        try:
            process_qc()
        except:
            print('error processing qc')


def process_qc(path='dflt', fails=False):
    summarize_qc(path=path)
    cleanup_qc(path=path)
    distribute_qc(path=path)
    resize_qc(path=path)
    if fails is True:
        combine_fails(path=path)


def summarize_qc(path='dflt'):
    """ pulls together qc data for all qc folders in path, and outputs a concise summary as well
    as a comprehensive consolidation of all qc metrics for each plate/batch """
    if path is 'dflt':
        path = gt.dflt_outpath(fldr_name='newQC')

    names, plates = [], []
    # n/a process control rows, only 1st column has volues
    rows = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 33, 34,
            65, 66, 67, 68, 69, 70, 71, 76, 77, 78, 79, 80,
            81, 82, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96]

    # load better data headers from reference text file
    headers = ['Range REF', 'Range PosAmp', 'Span REF', 'Span PosAmp',
        'FlogP REF', 'FlogP PosAmp', 'IQR REF', 'IQR PosAmp', 'Slope REF',
        'Slope PosAmp', 'Range plate avg', 'Range plate %cv', 'Span plate avg',
        'Span plate %cv', 'FlogP plate avg', 'FlogP plate %cv', 'IQR plate avg',
        'IQR plate %cv', 'Slope plate avg', 'Slope plate %cv', 'Beads median count',
        'Beads plate %cv', 'wells total', 'wells pass', 'wells fail',
        'vehicle wells total', 'vehicle wells pass', 'vehicle wells fail', 'MAD = 0',
        '0.02>= MAD >0', '0.05>= MAD >0.02', '0.1>= MAD >0.05', '0.2>= MAD >0.1',
        'Ref Affy Pear', 'Ref Affy Sear', 'Samplevar vehicle', 'Samplevar vehicle std',
        'Samplevar test', 'Samplevar test std', 'Samplevar poscon', 'Samplevar poscon std',
        '1st percent - MFI', 'level1 - MFI', 'level2 - MFI', 'level3 - MFI', 'level4 - MFI',
        'level5 - MFI', 'level6 - MFI', 'level7 - MFI', 'level8 - MFI', 'level9 - MFI',
        'level10 - MFI', 'level1 - cv', 'level2 - cv', 'level3 - cv', 'level4 - cv',
        'level5 - cv', 'level6 - cv', 'level7 - cv', 'level8 - cv', 'level9 - cv',
        'level10 - cv', 'cell ID', 'cell distance', 'MFI dist 95-5', 'Count dp52 median',
        'Count dp52 well std', 'Count dp52 bead std', 'dp52 rMCF7 lvl1', 'dp52 rMCF7 lvl10',
        'dp52 pMCF7 lvl1', 'dp52 pMCF7 lvl10', 'dp52 plateavg lvl1', 'dp52 plateavg lvl10',
        'dp52 5th pcnt    lvl10', 'dp52 95th pcnt    lvl10', 'Count dp53 median', 'Count dp53 well std',
        'Count dp53 bead std', 'dp53 rMCF7 lvl1', 'dp53 rMCF7 lvl10', 'dp53 pMCF7 lvl1', 'dp53 pMCF7 lvl10',
        'dp53 plateavg lvl1', 'dp53 plateavg lvl10', 'dp53 5th pcnt    lvl10', 'dp53 95th pcnt    lvl10',
        'dp52 - weak regions', 'dp52 lowcount', 'dp53 - weak regions', 'dp53 lowcount', '52 serial',
        '52 det date', '52 det duration', '53 serial', '53 det date', '53 det duration']

    # use '**/*' with glob recursive=True to go through subfolders
    searcht = os.path.join(path, '**/*QC_metrics.txt')
    for i, file in enumerate(glob.iglob(searcht, recursive=True)):
        m = pd.read_csv(file, sep='\t')
        m['index'] = headers
        m = m.set_index('index')
        # drop blank columns (all but batches w/ data)
        m = m[m.columns[~m.columns.str.contains('Unnamed')]]
        # fill across n/a values of process controls to other batches
        if len(m.columns) > 1:
            for col in range(1,len(m.columns)):
                m.iloc[rows,col] = m.iloc[rows,0]
        # build up list of all batches in folder
        name = m.columns.values[0]
        try:
            nf = file.replace('QC_metrics', 'highcv')
            cv = pd.read_csv(nf, sep='\t')
            cv50 = round(cv.iloc[49, 2], 2)
            cv500 = round(cv.iloc[499, 2], 2)
            nf = file.replace('QC_metrics', 'periodicity_dot')
            pdcty = pd.read_csv(nf, sep='\t')
            pdcty50 = round(pdcty.iloc[49, 2], 2)
            pdcty150 = round(pdcty.iloc[149, 2], 2)
            m.loc['cv50'] = [cv50] * len(m.columns)
            m.loc['cv500'] = [cv500] * len(m.columns)
            m.loc['pdcty50'] = [pdcty50] * len(m.columns)
            m.loc['pdcty150'] = [pdcty150] * len(m.columns)
        except FileNotFoundError:
            print('no highcv/periodicity.txt found')
            m.loc['cv50'] = ['na'] * len(m.columns)
            m.loc['cv500'] = ['na'] * len(m.columns)
            m.loc['pdcty50'] = ['na'] * len(m.columns)
            m.loc['pdcty150'] = ['na'] * len(m.columns)
        plates.append(m)

    metrics = pd.concat(plates, axis=1)

    # merge all metrics into one dataframe
    metrics = metrics.T
    metrics.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    # define output file name
    opath = os.path.join(path, os.path.basename(path) + '_qcmetrics.txt')
    metrics.rename(columns = {'Unnamed: 0': 'plate'}, inplace=True)
    # if all single batch plates, remove batch identifier
    if any('-B' in x for x in metrics.index.values):
        pass
    else:
        metrics.index = metrics.index.str.replace('-A', '')
    metrics.sort_index(inplace=True)
    metrics.to_csv(opath, sep='\t')

    # pipeline output headers
    orighdr = ['level1 - MFI', 'Ref Affy Pear', 'FlogP REF', 'FlogP plate avg',
        'cell distance', 'cv50', 'pdcty50', 'Samplevar vehicle', 'Samplevar test',
        'wells fail', 'vehicle wells pass']
    # the headers i want
    newhdr = ['Plate L1', 'Ref Acorr', 'Ref Flogp', 'Plate Flogp', 'cell dist', 'cv50',
        'pdcty50', 'Veh var', 'Test var', 'failed wells', 'vehicle pass']
    # the new headers I derive
    addnlhdr = ['Pos L10', 'Ref L10', 'Plate L10', 'Plate L10 var']
    #create subset
    try:
        subset = metrics[orighdr]
        subset.columns = newhdr
    except KeyError:
        orighdr = ['level1 - MFI', 'Ref Affy Pear', 'FlogP REF', 'FlogP plate avg',
        'cell distance', 'Samplevar vehicle', 'Samplevar test',
        'wells fail', 'vehicle wells pass']
        newhdr = ['Plate L1', 'Ref Acorr', 'Ref Flogp', 'Plate Flogp', 'cell dist',
        'Veh var', 'Test var', 'failed wells', 'vehicle pass']
        subset = metrics[orighdr]
        subset.columns = newhdr
    # derive the additional summary numbers from existing
    # apply pd.to_numeric because it's being flaky and not doing things right
    lvl10 = metrics[['dp52 plateavg lvl10', 'dp53 plateavg lvl10']].apply(lambda x: pd.to_numeric(x, errors='ignore'))
    s10 = lvl10.mean(axis=1)
    r10 = metrics[['dp52 rMCF7 lvl10', 'dp53 rMCF7 lvl10']].apply(lambda x: pd.to_numeric(x, errors='ignore'))
    r10 = r10.mean(axis=1)
    p10 = metrics[['dp52 pMCF7 lvl10', 'dp53 pMCF7 lvl10']].apply(lambda x: pd.to_numeric(x, errors='ignore'))
    p10 = p10.mean(axis=1)
    lvl10cv = metrics['level10 - cv'].apply(lambda x: pd.to_numeric(x, errors='ignore'))
    p10var = lvl10cv / s10
    p10var.dropna(inplace=True)
    p10var = p10var.map('{:,.2f}'.format)
    #add into the complete summary file
    addnl = pd.concat([p10, r10, s10, p10var], axis=1)
    addnl.columns = addnlhdr
    full = pd.concat([addnl, subset], axis=1)
    opath = opath.replace('metrics', 'summary')
    full.sort_index(inplace=True)
    full.to_csv(opath, sep='\t')
    # additionall combine qc fails into single file
    combine_fails(path=path)
    try:
        run_ref_pca(path)
    except ValueError:
        pass


def cleanup_qc(path='dflt'):
    if path is 'dflt':
        path = gt.dflt_outpath(fldr_name='QCprocessing')

    #for folder in folders:
    shutil.rmtree(path)


def distribute_qc(path = 'dflt'):
    if path is 'dflt':
        inpath = gt.dflt_outpath(fldr_name='newQC')
        outpath = gt.dflt_outpath(fldr_name='QCprocessing')

    folders = ['calibs', 'flogps', 'escore', 'cellid-nolabel', 'cellid-label', 'euclidean']
    folders = [os.path.join(outpath, x) for x in folders]
    srch_terms = ['finalqc/*calibplot', 'finalqc/*FLOGP', 'escore_summary*/', '*cell_line/*cellid_nolabel/*-*cellid_circle',
                  '*cell_line/*-*cellid_circle', '-*euclidean']

    for term, fold in zip(srch_terms, folders):
        try:
            os.makedirs(fold)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
        srch = '*'.join(['', term, '' ]) + '.png'
        for file in gt.globit(inpath, srch):
            shutil.copy(file, fold)


def resize_qc(path='dflt'):
    """ run through folder contents and subdirs and conduct appropriate image resizing"""
    if path is 'dflt':
        path = gt.dflt_outpath(fldr_name='QCprocessing')
    fl = glob.glob(path + '/**/*.png', recursive=True)
    for f in fl:
        if 'calibplot' in f:
            calplot(f, f)
        elif '_es1' in f:
            escoresum(f, f)
        elif 'FLOGP' in f:
            flogp(f, f)
        elif 'euclidean' in f:
            euclidean(f, f)
        elif 'cellid' in f:
            cellid(f, f)


# sub processes below

def resize_batch_calplots(path='dflt', outpath='dflt', plot_type='clean'):
    if path is 'dflt':
        path = '/Users/WRB/Desktop/newQC/'
    if outpath is 'dflt':
        outpath = '/Users/WRB/Desktop/newfigs/'
        if not os.path.exists(outpath):
            os.mkdir(outpath)

    plate_names = [gt.get_shn(x) for x in os.listdir(path) if x[0] is not '.']
    plate_names = set([x for x in plate_names if len(x)==6])

    for pn in plate_names:
        calfold = path + pn + '_qc/' + pn + '_batch_figs/' + pn + '_calibplot/'
        for batch in ['A', 'B', 'C', 'D']:
            loc = calfold + batch + '/'
            fn = pn + '_' + plot_type + '_calibplot_dp'
            for bset in ['52','53']:
                file = fn + bset + '.png'
                filepath = (outpath + file).replace(plot_type, batch)
                calplot(loc + file, filepath, crop=True)


def barviews(path):
    fl = glob.glob(path + '/**/*barview*.png', recursive=True)
    for f in fl:
        im = Image.open(f)
        new_im = im.resize((69, 689), Image.ANTIALIAS)
        new_im.save(f, dpi=(220,220))


def calplot(file, outpath, crop=False):
    im = Image.open(file)
    if crop is True:
        im_crop = im.crop((307, 133, 2175, 1605))
        #im_crop = im.crop((618, 252, 4362, 3217))
    else:
        im_crop = im
    im_resized = im_crop.resize((781,616), Image.ANTIALIAS)
    if outpath.endswith('.png'):
        pass
    else:
        outpath += '.png'
    im_resized.save(outpath, dpi=(220, 220))


def cellid(file, outpath):
    im = Image.open(file)
    #im_resized = im.resize((484,484), Image.ANTIALIAS)
    im_resized = im.resize((440, 440), Image.ANTIALIAS)
    if outpath.endswith('.png'):
        pass
    else:
        outpath += '.png'
    im_resized.save(outpath, dpi=(220, 220))


def flogp(file, outpath):
    im = Image.open(file)
    im_resized = im.resize((896,569), Image.ANTIALIAS)
    if outpath.endswith('.png'):
        pass
    else:
        outpath += '.png'
    im_resized.save(outpath, dpi=(220, 220))


def escoresum(file, outpath):
    im = Image.open(file)
    new_width = int((880 * im.size[0]) / im.size[1])
    im_resized = im.resize((new_width,880), Image.ANTIALIAS)
    if outpath.endswith('.png'):
        pass
    else:
        outpath += '.png'
    im_resized.save(outpath, dpi=(220, 220))


def euclidean(file, outpath):
    im = Image.open(file)
    im_crop = im.crop((160, 0, 749, 476))
    im_resized = im_crop.resize((722,581), Image.ANTIALIAS)
    if outpath.endswith('.png'):
        pass
    else:
        outpath += '.png'
    im_resized.save(outpath, dpi=(220, 220))


def sub_check_lines(g):
    """ with an (inf) gct file path will check number of rows corresponds with gct header """
    i = 0
    with open(g.file, 'r') as f:
        for i, l in enumerate(f):
            i += 1
    lines = g.drows + g.hrows + 3
    if lines == i:
        return True
    else:
        print('line error i={} lines={}'.format(i, lines))
        return False


def sub_check_columns(g):
    """ with an (inf) gct file path will check number of columns corresponds with gct header """
    with open(g.file, 'r') as f:
        maxlen = 0
        for l in f:
            length = len(l.split('\t'))
            if length > maxlen:
                maxlen = length
    cols = g.dcols + g.hcols + 1
    if cols == maxlen:
        return True
    else:
        print('error maxlen={} cols={}'.format(maxlen, cols))
        return False


def sub_check_failed(g, txt):
    """ provided (inf) gct file and final map txt file checks if the failed wells line up
     currently just counts not well identities """
    fails = 0; empty = 0
    try:
        with open(txt, 'r') as f:
            for l in f:
                if 'failed' in l:
                    fails += 1
                elif 'empty' in l:
                    empty += 1
    except:
        try:
            m = pd.read_excel(txt)
            fails = len(m[m['type']=='failed'])
            empty = len(m[m['type'] == 'empty'])
        except:
            pass

    if g.dcols == 384 - fails - empty:
        return fails, True
    else:
        print('fail error fails={} cols={}'.format(fails, g.dcols))
        return False


def sub_check_decimal(g):
    """ check number of data decimal places """
    with open(g.file, 'r') as file:
        for i in range(5):
            file.readline()
        line = file.readline()
        val = line.split('\t')[3]
        dplaces = len(val.split('.')[1])
    return dplaces


def main():
    pass


if __name__ == '__main__':
    main()

from PIL import Image
import pt, os, glob, gt, shutil, errno, boto3, subprocess
import pandas as pd


def combine_fails(path='dflt'):
    if path == 'dflt':
        path = gt.dflt_outpath(fldr_name='newQC')
    fl = gt.globit(path, '*QC_fail*')
    files = ' '.join(fl)
    cmd_str = 'cat ' + files + ' > ' + os.path.join(path, 'QC_fail.txt')
    subprocess.run(cmd_str, shell=True)


def dl_data(sc='q', src='dflt', dest='dflt', clr=False, search=None, excl=None, ext=None, p=False):

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
        s3c.download_file('genometry', 'PCA_analysis/PCA.pptx', os.path.join(dest, 'PCA.pptx'))
        coords = s3c.list_objects(Bucket='genometry', Prefix='PCA_analysis/')
        coord = sorted([x['Key'] for x in coords['Contents']])[-1]
        s3c.download_file('genometry', coord, os.path.join(dest, 'PCA_coords.txt'))
    else:
        src = 's3://genometry/' + src

    search_args = []

    if search is None:
        if 'q' in sc:
            search_args.append(('*_qc/*', ''))
        if 'g' in sc:
            search_args.append((['*_fullqnorm_*', '*_QNORM_sorted*'], ''))
        if 'z' in sc:
            search_args.append(('*_ZSVCQNORM_*', ''))
        if 'e' in sc:
            search_args.append(('*_escore/*', '*.gct'))
    else:
        if '*' not in search:
            search = '*' + search + '*'
        search_args = [(search, excl)]

    for search, excl in search_args:
        cmd_str = f'aws s3 cp --recursive {src} {dest} --exclude "*"'
        if isinstance(search, str):
            search = [search]
        for st in search:
            cmd_str += f' --include {st}'
        if excl != '':
            cmd_str += f' --exclude {excl}'

        print(cmd_str)

        subprocess.run(cmd_str, shell=True)

    if ext is not None:
        fl = gt.globit(dest, f'*{ext}')
        for f in fl:
            file_dest = os.path.join(dest, os.path.basename(f))
            shutil.move(f, file_dest)
        subdirs = [x[0] for x in os.walk(dest)][1:]
        for sd in subdirs:
            os.rmdir(sd)

    if p is True:
        process_qc()


def process_qc(path='dflt'):
    summarize_qc(path=path)
    combine_fails(path=path)
    cleanup_qc(path=path)
    distribute_qc(path=path)
    resize_qc(path=path)


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
    lvl10 = metrics[['dp52 plateavg lvl10', 'dp53 plateavg lvl10']].apply(pd.to_numeric)
    s10 = lvl10.mean(axis=1)
    r10 = metrics[['dp52 rMCF7 lvl10', 'dp53 rMCF7 lvl10']].apply(pd.to_numeric)
    r10 = r10.mean(axis=1)
    p10 = metrics[['dp52 pMCF7 lvl10', 'dp53 pMCF7 lvl10']].apply(pd.to_numeric)
    p10 = p10.mean(axis=1)
    lvl10cv = metrics['level10 - cv'].apply(pd.to_numeric)
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


def cleanup_qc(path='dflt'):
    if path is 'dflt':
        path = gt.dflt_outpath(fldr_name='QCprocessing')

    #for folder in folders:
    shutil.rmtree(path)


def distribute_qc(path = 'dflt'):
    if path is 'dflt':
        inpath = gt.dflt_outpath(fldr_name='newQC')
        outpath = gt.dflt_outpath(fldr_name='QCprocessing')

    folders = ['calibs', 'flogps', 'escore', 'cellid', 'euclidean']
    folders = [os.path.join(outpath, x) for x in folders]
    srch_terms = ['finalqc/*calibplot', 'finalqc/*FLOGP', 'escore_summary*/', '-*cellid_circle', '-*euclidean']

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


def resize_batch_calplots(path='dflt', outpath='dflt', plot_type='clean'):
    if path is 'dflt':
        path = '/Users/WRB/Desktop/newQC/'
    if outpath is 'dflt':
        outpath = '/Users/WRB/Desktop/newfigs/'
        if not os.path.exists(outpath):
            os.mkdir(outpath)

    plate_names = [pt.get_shn(x) for x in os.listdir(path) if x[0] is not '.']
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
    im_resized = im.resize((484,484), Image.ANTIALIAS)
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


def main():
    pass


if __name__ == '__main__':
    main()
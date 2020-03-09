#!/usr/bin/python


import collections as cll
import pandas as pd
import gct, os, sys, gt, pa
import matplotlib.pyplot as plt

def enrich_matrix(df, h,labels='dflt',type='scaled'):
    """ from a passed dataframe, calculate pairwise enrichment
    auto thresholds of features above zscore limited by number """
    zsthresh = 2
    numthresh = 50
    en_mtrx = pd.DataFrame()
    labels = gt.check_dfltarg(labels, h.label)
    for sample, label in zip(df.columns, h.label.values):
        # determine thresholds
        (up, dn) = sigs.get_sig(df[sample], zsthresh, numlim=numthresh)

        enres = sigs.bulk_test_enrich((up, dn), df, h)

        en_mtrx = pd.concat([en_mtrx, enres[type]], axis=1)

    en_mtrx.columns = h.label
    en_mtrx.set_index(h.label, inplace=True)
    return en_mtrx


def signature_overlap_matrix(df, labels, thresh=2):
    """ given a sorted dataframe, derive signatures from each instance and plot separate up and down
    overlap matrices of signature genes """
    uplist, downlist = [], []
    for sample in df.columns:
        up, dn = get_sig(df[sample], thresh)
        uplist.append(up)
        downlist.append(dn)
    up_matrix = gt.overlap_matrix(uplist, labels)
    dn_matrix = gt.overlap_matrix(downlist, labels)
    return up_matrix, dn_matrix


def load_sig(fpaths):
    """ load tuple sig of two files, either as one string or list of two paths """
    if isinstance(fpaths, str):
        pathlist = gt.splitpaths(fpaths, ext='.grp')
        up_path = [x for x in pathlist if '_up' in x][0]
        dn_path = [x for x in pathlist if '_dn' in x][0]
    else:
        try:
            # print('signature load {} paths.'.format(len(fpaths)))
            up_path = [x for x in fpaths if '_up' in x][0]
            dn_path = [x for x in fpaths if '_dn' in x][0]
        except:
            print('error with load sig input')
            pathlist = ['foo','foo']

    with open(up_path, 'r') as file:
        up = file.readlines()
        up = [x.strip() for x in up]

    with open(dn_path, 'r') as file2:
        dn = file2.readlines()
        dn = [x.strip() for x in dn]


    return (up, dn)


def convert_to_symbols(probelist):
    """ convert affx probe ids to gene symbols. currently hardcoded the location of ref file """
    dict_loc = '/Users/WRB/Dropbox/Areas of Focus/_Genometry/Oligos/Genometry L1000 Genes (25 March 2015).xlsx'
    gd = pd.read_excel(dict_loc, skiprows=4, index_col=0)
    gd = gd.iloc[:, 0:6]
    genelist = []
    for p in probelist:
        try:
            genelist.append(gd.loc[p, 'Gene Symbol'])
        except:
            print('problem with lookup ', p)
    return genelist


def sig_survey(inst, cust=None):
    """ run to quickly get an overview of signature sizes with different zs cutoffs """
    sigdict = cll.defaultdict(list)
    threshs = [1.5, 2, 3, 5, 8, 10]

    if isinstance(inst, pd.DataFrame):
        print('getting consensus for passed dframe')
        inst = pa.consensus(inst)

    if cust:
        if gt.isnum(cust):
            threshs.append(float(cust))
        elif len(cust) > 1:
            for num in cust:
                threshs.append(float(num))
    try:
        inst = inst[inst.columns[0]]
    except AttributeError:
        pass

    for t in threshs:
        sigdict[t].append(len(inst[inst >= t]))
        sigdict[t].append(len(inst[inst <= -t]))
    result = ''
    for k, v in sigdict.items():
        result += ('{}:({}/{})  '.format(k, v[0], v[1]))
    return result


def get_sig(inst, t, td=None, numlim=None):
    """ put in pd.series of zscore values with zscore threshold, returns tuple of (up,down) lists of genes
    if different thresholds are desired, t will refer to up threshold and td is the down threshold """

    if isinstance(inst, pd.DataFrame):
        print('getting consensus for passed dframe')
        inst = pa.consensus(inst)

    inst.sort_values(ascending=False, inplace=True)

    up = inst[inst >= t].index.values

    if td is None:
        dn = inst[inst <= -t].index.values
    else:
        dn = inst[inst <= -td].index.values

    if numlim is not None:
        if len(up) > numlim:
            up = up[:numlim]
        if len(dn) > numlim:
            dn = dn[:numlim]

    return (up, dn)


def save_sig(sig, gs=False, name='dflt', path='dflt'):
    """ pass a signature tuple of up/down and name/path to save those lists """
    if name is 'dflt':
        name = 'mysig'
        try:
            name = name.replace(':','-')
        except AttributeError:
            pass
    if path is 'dflt':
        path = gt.dflt_outpath(fldr_name='foo')
    up, dn = sig[0], sig[1]
    if gs is True:
        up = convert_to_symbols(up)
        dn = convert_to_symbols(dn)
    filename = name + '_up.grp'
    gt.savelist(up, os.path.join(path, filename))
    filename = name + '_dn.grp'
    gt.savelist(dn, os.path.join(path, filename))


def get_and_save_sig(inst, t, gs=False, name='dflt', path='dflt'):
    """ automatically save generated signature with default or provided file name
    and destination folder. gs flag saves things in terms of gene symbols """
    if name is 'dflt':
        name = inst.name
        try:
            name = name.replace(':','-')
        except AttributeError:
            pass
    if path is 'dflt':
        path = gt.dflt_outpath(fldr_name='foo')
    up, dn = get_sig(inst, t)
    if gs is True:
        up = convert_to_symbols(up)
        dn = convert_to_symbols(dn)
    filename = name + '_up.grp'
    gt.savelist(up, os.path.join(path, filename))
    filename = name + '_dn.grp'
    gt.savelist(dn, os.path.join(path, filename))
    return (up, dn)


def bulk_test_enrich(sig, df, h, outpath=False):
    """ pass in a dataframe and a signature (tuple of up/down), and optionally map/header
    information to include in the returned enrichment scores"""
    up, dn = sig[0], sig[1]
    escore = {}
    for c in df.columns:
        escore[c] = test_enrichment(df[c], up, dn)
    edf = pd.DataFrame(escore)
    edf = edf.T
    # create local scaled enrichment
    pmax = edf['absolute'].max()
    pmin = edf['absolute'].min()
    edf['scaled'] = edf['absolute'].apply(lambda x: x/pmax if x > 0 else (x/(-1*pmin)) if x < 0 else 0)
    edf['scaled'] = edf['scaled'].apply(lambda x: float('{:.3f}'.format(x)))
    edf = edf[['scaled', 'absolute', 'up', 'dn']]
    # optionally merge results with sample header obj
    if h is not None:
        edf = pd.merge(h, edf, left_index=True, right_index=True, how='inner')
    edf.sort_values('scaled', ascending=False, inplace=True)
    if outpath is not False:
        if outpath == 'dflt':
            outpath = gt.dflt_outpath(fn=df.name + '_enrichment.xlsx')
        edf.to_excel(outpath)
    return edf


def make_barview(edf, argdict, ax=None, label=False, height=2):
    """ given an edf header file with enrichment info, and an argument dictionary for highlighted
    intances within the barview. passed label will be on left hand side

    the height of the highlighted instances can be sensitive to being swaamped out and invisible

    """

    my_ax = ax

    ax = format_barview_plot(edf, ax=my_ax)

    if label is True:
        ax.set_ylabel(list(argdict.values())[0], labelpad=0.0)
    elif label is not False:
        ax.set_ylabel(label, labelpad=0.0)

    edf.sort_values(['scaled', 'up'], ascending=False, inplace=True)

    edf['plot_pos'] = list(range(len(edf)+2, 2, -1))
    edf['rank'] = list(range(1, len(edf)+1, 1))

    pos = edf[edf['scaled'] > 0]
    null = edf[edf['scaled'] == 0]
    neg = edf[edf['scaled'] < 0]
    selected = gt.hsub(edf, argdict)

    sslist = [pos, null, neg, selected]
    clist = ['lime', 'lightgrey','red', 'black']

    # height of the highlight bar can be sensitive
    for subset, color in zip(sslist, clist):
        ax.barh(subset['plot_pos'], [1*len(subset)], color=color, align='center', height=height)

    if my_ax is None:
        plt.tight_layout()


def make_barview_range(edf, argdict, across='dose', label=False, outpath=False):
    """ with enrichment score results, plot barviews across the range of conditions, default dose """

    cond_range = sorted(gt.hsub(edf, argdict)[across].unique())

    print(argdict.values())
    mytitle = ' '.join(argdict.values())

    fig, axarr = plt.subplots(1, len(cond_range), sharey='row')

    for i,cond in enumerate(cond_range):
        my_ax = axarr[i]
        new_argdict = argdict
        if across is not None:
            new_argdict[across] = cond
        if label is True:
            make_barview(edf, new_argdict, ax=my_ax, label=cond)
        else:
            make_barview(edf, new_argdict, ax=my_ax)

    #fig.subplots_adjust(hspace=0.5)

    plt.suptitle(mytitle)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    if outpath is True:
        outpath = gt.dflt_outpath(fldr_name='foo')
        myoutpath = os.path.join(outpath, mytitle +'_enrich.png')
        plt.savefig(myoutpath)
        plt.close()


def format_barview_plot(edf, ax=None):
    """ plot a barview, with edf passed in or just...? """
    if ax is None:
        fig = plt.figure(figsize=[.3 , 3])
        ax = fig.add_subplot(111)
    ax.set_facecolor("white")
    ax.set_ylim([0, len(edf)+3])
    ax.set_xlim([0, 1])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.tick_params(axis='x', bottom='off', top='off', labelbottom='off')
    ax.tick_params(axis='y', right='off', left='off')
    ax.set_yticklabels([])
    ax.yaxis.set_label_position('left')
    return ax


def test_enrichment(inst, up, dn, v=False):
    """ signature enrichment implemented from the cmap help how to documents
    calculates and return the up, down and absolute scores for the given up and dn
    tag files and the pd.Series object instance
    v is for Verbose """
    try:
        inst = inst.sort_values(ascending=False)
    except:
        inst = inst[inst.columns[0]].sort_values(ascending=False)
    igenes = list(inst.index.values)
    escore = {}

    for nm, sig in zip(['up', 'dn'], [up, dn]):
        ascore, bscore = 0, 0

        my_sig = sorted(sig, key=lambda x: igenes.index(x))

        if v is True:
            print([igenes.index(x) for x in my_sig])

        for i, g in enumerate(my_sig, 1):
            # test to make sure gene is in index, quit otherwise
            try:
                igenes.index(g)
            except ValueError:
                sys.exit('gene '+ g + ' not in instance index len=' + str(len(igenes)))

            # print('sig gene #' + str(i) + ' - ' + g + ' is in position ' + str(igenes.index(g)))

            aval = (i / len(my_sig)) - (igenes.index(g) / len(igenes))

            bval = (igenes.index(g)/len(igenes)) - ((i - 1) / len(my_sig))

            if aval > ascore:
                ascore = aval
            if bval > bscore:
                bscore = bval

        if v is True:
            print(f'ascore is {ascore:.3}, and bscore is {bscore:.3}')

        if ascore >= bscore:
            subscore = ascore
        elif ascore < bscore:
            subscore = -1 * bscore
        escore[nm] = round(subscore,3)
    if escore['up'] > 0 and escore['dn'] > 0:
        escore['absolute'] = 0
    elif escore['up'] < 0 and escore['dn'] < 0:
        escore['absolute'] = 0
    else:
        escore['absolute'] = escore['up'] - escore['dn']
    return escore


def get_s2n_genes(g, c1, c2):
    d, h = gct.extractgct(g)
    c1w = gt.hsub(h, {'well':c1}).index
    c2w = gt.hsub(h, {'well':c2}).index
    d1 = d[c1w]
    d2 = d[c2w]
    res = sig_to_noise(d1, d2)
    return res


def sig_to_noise(coh1, coh2, drop=False):
    """ signal to noise calculation for looking at features differentially
    expressed between two cohorts of samples, where cohorts are dataframes """
    results = pd.DataFrame(index=coh1.index)
    for i, c in enumerate([coh1, coh2]):
        results['mean ' + str(i)] = c.mean(axis=1)
        results['std ' + str(i)] = c.std(axis=1)
    print(coh1.head())

    results['s2n'] = (results['mean 0'] - results['mean 1'])/(results['std 0'] + results['std 1'])
    results['abs_s2n'] = abs(results['s2n'])
    results.sort_values('abs_s2n', ascending=False, inplace=True)
    if drop is True:
        results.dropna(inplace=True)
    return results


def survey_granularity(df, c=None):
    vals = df.apply(lambda x: len(x.unique()), axis=1)
    vals = pd.DataFrame(vals, columns=['num_vals'])
    vals['med'] = df.apply(lambda x: x.median(), axis=1).values
    vals.sort_values('num_vals', inplace=True)
    print(len(vals[vals['num_vals'] <= 30]))
    if c != None:
        sv = vals[vals['num_vals'] <= 30]
        c.update(sv.index.values)
    else:
        return vals


def run_granularity(path):
    flist = gt.get_flist(path, '.gct')
    c = cll.Counter()
    for f in flist:
        d, h = gct.extractgct(f)
        survey_granularity(d, c)
    c = pd.Series(c, name='count')
    c.sort_values(ascending=False, inplace=True)
    c = c[c > 1]
    c.to_excel(os.path.join(path, 'counter.xlsx'))





def main():
    pass


if __name__ == '__main__':
    main()



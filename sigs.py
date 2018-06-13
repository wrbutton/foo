#!/usr/bin/python


import collections as cll
import pandas as pd
import pt, gct, os, sys, gt


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
    threshs = [2, 3, 5, 8, 10]
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


def get_sig(inst, t):
    """ put in pd.series of zscore values with zscore threshold, returns up,down lists"""
    try:
        inst = inst[inst.columns[0]]
    except AttributeError:
        pass
    up = inst[inst >= t].index.values
    dn = inst[inst <= -t].index.values
    return up, dn


def save_sig(inst, t, gs=False, name='dflt', path='dflt'):
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


def bulk_test_enrich(df, up, dn, h=None, outpath=False):
    """ pass in a dataframe and a pair of signature lists, and optionally map/header
    information to include in the returned enrichment scores"""
    escore = {}
    for c in df.columns:
        escore[c] = test_enrichment(df[c], up, dn)
    edf = pd.DataFrame(escore)
    edf = edf.T
    # create local scaled enrichment
    pmax = edf['abslt'].max()
    pmin = edf['abslt'].min()
    edf['scaled'] = edf['abslt'].apply(lambda x: x/pmax if x > 0 else (x/(-1*pmin)) if x < 0 else 0)
    edf['scaled'] = edf['scaled'].apply(lambda x: '{:.3f}'.format(x))
    edf = edf[['scaled', 'abslt', 'up', 'dn']]
    # optionally merge results with sample header obj
    if h is not None:
        edf = pd.merge(h, edf, left_index=True, right_index=True, how='inner')
    if outpath is not False:
        if outpath == 'dflt':
            outpath = gt.dflt_outpath(fn=df.name + '_enrichment.xlsx')
        edf.to_excel(outpath)
    return edf


def make_barview(df, up, dn, argdict):
    """ given an edf header file with enrichment info, and an argument dictionary """


def test_enrichment(inst, up, dn):
    """ signature enrichment implemented from the cmap help how to documents
    calculates and return the up, down and absolute scores for the given up and dn
    tag files and the pd.Series object instance """
    igenes = list(inst.sort_values(ascending=False).index.values)
    escore = {}
    for nm, sig in zip(['up', 'dn'], [up, dn]):
        ascore, bscore = 0, 0
        for i, g in enumerate(sig, 1):
            try:
                aval = (i / len(sig)) - (igenes.index(g)/len(igenes))
            except ValueError:
                sys.exit('gene '+ g+ ' not in instance index len=' + str(len(igenes)))
            # print('sig gene #' + str(i) + ' - ' + g + ' is in position ' + str(igenes.index(g)))
            bval = (igenes.index(g)/len(igenes)) - ((i - 1) / len(sig))
            if aval > ascore:
                ascore = aval
            if bval > bscore:
                bscore = bval
        if ascore >= bscore:
            subscore = ascore
        elif ascore < bscore:
            subscore = -1 * bscore
        escore[nm] = round(subscore,3)
    if escore['up'] > 0 and escore['dn'] > 0:
        escore['abslt'] = 0
    elif escore['up'] < 0 and escore['dn'] < 0:
        escore['abslt'] = 0
    else:
        escore['abslt'] = escore['up'] - escore['dn']
    return escore


def get_s2n_genes(g, c1, c2):
    d, h = gct.extractgct(g)
    #c1w = h[h['name']==c1]
    #c2w = h[h['name']==c2]
    c1w = gt.hsub(h, {'well':c1}).index
    c2w = gt.hsub(h, {'well':c2}).index
    d1 = d[c1w]
    d2 = d[c2w]
    res = sig_to_noise(d1, d2)
    return res
    #genes = res[res['abs_s2n'] > t]


def sig_to_noise(coh1, coh2, drop=False):
    """ signal to noise calculation for looking at features differentially
    expressed between two cohorts of samples """
    for c in [coh1, coh2]:
        c['mean'] = c.mean(axis=1)
        c['std'] = c.std(axis=1)
    print(coh1.head())
    results = pd.DataFrame(index=coh1.index)
    results['s2n'] = (coh1['mean'] - coh2['mean'])/(coh1['std'] + coh2['std'])
    results['abs_s2n'] = abs(results['s2n'])
    results.sort_values('s2n', inplace=True)
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
    flist = pt.get_flist(path, '.gct')
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



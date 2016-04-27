#!/usr/bin/env python

from functools import partial
import logging
import multiprocessing as mp
import re
import time

import pandas as pd
import numpy as np
from sqlalchemy import MetaData, create_engine
from scipy.stats import pearsonr
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import ttest_ind

import leip

lg = logging.getLogger(__name__)

@leip.arg('-s', '--set', action='append', nargs=2, help='set to value')
@leip.arg('-r', '--keep', default='.*', help='regex - retain these samples')
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('output_file')
@leip.command
def generate_template(app, args):
    engine = create_engine(args.datafile)
    lg.info('load from: %s', args.datafile)
    c = pd.read_sql('select * from normcounts limit 5',
                        engine, index_col='index')
    samples = c.columns
    if args.keep:
        rekeep = re.compile(args.keep)
        samples = list(filter(rekeep.search, samples))
        
    template = pd.DataFrame([0] * len(samples), index=samples)
    
    for setpat, setval in args.set:
        setval = float(setval)
        reset = re.compile(setpat)
        setlist = filter(reset.search, samples)
        template.loc[setlist] = setval
    template.to_csv(args.output_file, sep="\t", header=None)


def ttest(row, setA, setB):
    row = row[1]
    return ttest_ind(row[setA], row[setB])

def apply_pearson(r, tmpl):
    name, r = r
    return pearsonr(r, tmpl)

def apply_euclid(r, tmpl):
    return 1 - (2 * (np.sum(np.abs(r - tmpl)) / len(r)))

@leip.arg('-f', '--filter', action='append', nargs=3)
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('-j', '--threads', default=8, type=int)
@leip.arg('-m', '--method', help='pearson, euclid', default='pearson')
@leip.flag('--debug')
@leip.flag('-x', '--exponent', help='convert normalized count values '
           'from log space to count space')
@leip.arg('template')
@leip.command
def correlate(app, args):
    start_time = time.time()
    
    DEBUG = args.debug
    
    engine = create_engine(args.datafile)
    assert args.method in ['pearson', 'euclid']

    lg.info('loading template: %s', args.template)
    tmpl = pd.read_csv(args.template, sep="\t", index_col=0, header=None)[1]

    if DEBUG:
        lg.info('load fraction & stats data from: %s', args.datafile)
        d = pd.read_sql('select * from fraction_inf LIMIT 2500', engine, index_col='index')
        stats = pd.read_sql('select * from junction_stats LIMIT 2500', engine,
                            index_col='index')
    else:
        lg.info('load fraction & stats data from: %s', args.datafile)
        d = pd.read_sql('select * from fraction_inf', engine, index_col='index')
        stats = pd.read_sql('select * from junction_stats', engine,
                            index_col='index')

    #check if we're in binary mode - ie there are only two possibilities in the template
    binary = True if len(tmpl.value_counts()) == 2 else False
        
    if binary:
        lg.info("Binary mode on")

        lg.info("calc few fraction stats")
        group_values = list(tmpl.value_counts().index.to_series().sort_values(ascending=False))
        groupA, groupB = group_values
        sampleA = tmpl[tmpl == groupA].index
        sampleB = tmpl[tmpl == groupB].index
        meanA = d[sampleA].mean(1)
        meanB = d[sampleB].mean(1)
        fracdiff = meanA - meanB
        
        #calculate fractional differences
        lg.info("Calculate fraction difference")
        lg.info("-- min: %.3f", fracdiff.min())
        lg.info("-- max: %.3f", fracdiff.max())
        
        # ttest on fractions
        ttest2 = partial(ttest, setA=sampleA, setB=sampleB)
        with mp.Pool(args.threads) as P:
            rv = P.map(ttest2, d.iterrows())
            
        rv = pd.Series(rv)
        rv.index = d.index
        frac_tt = rv.str.get(0)
        frac_tp = rv.str.get(1).fillna(1)
        frac_tpadj = multipletests(frac_tp, method='fdr_bh')[1]
        frac_tpadj = pd.Series(frac_tpadj, index=d.index)
                
        lg.info("loading normalized counts")
        if DEBUG:
            norm = pd.read_sql('select * from normcounts LIMIT 5000', engine, index_col='index')
        else:
            norm = pd.read_sql('select * from normcounts', engine, index_col='index')
        lg.info("-- loaded %d records", len(norm))
        norm = 2 ** norm if args.exponent else norm
        assert norm.min().min() >= 0
        lg.info("-- start calculating LFC")
        lfc = np.log2(norm[sampleA].mean(1) / norm[sampleB].mean(1))
        lg.info("-- min lfc: %.3f", lfc.min())
        lg.info("-- max lfc: %.3f", lfc.max())

        # ttest on junctions
        with mp.Pool(args.threads) as P:
            rv = P.map(ttest2, norm.iterrows())
            
        rv = pd.Series(rv)
        rv.index = norm.index
        norm_tt = rv.str.get(0)
        norm_tp = rv.str.get(1)
        norm_tpadj = pd.Series(multipletests(norm_tp, method='fdr_bh')[1], index=norm.index)
        
    lg.info("raw data shape: %s", d.shape)
    lg.info("remove samples which are not in the template")
    d = d[tmpl.index]
    lg.info("resulting data shape: %s", d.shape)

    if args.filter:
        lg.info("start filtering %d junction records", len(stats))
    for fcol, fop, fcutoff in args.filter:
        sc = stats[fcol]
        flt = getattr(sc, fop)(float(fcutoff))
        stats = stats[flt]
        lg.info("after %s %s %s, %d records left", fcol, fop, fcutoff, len(stats))


    lg.info('template value counts:' + " ".join(str(tmpl.value_counts()).split()))
    
    lg.info("start calculating correlation, method: %s", args.method)

    if args.method == 'pearson':
        apply_pearson_2 = partial(apply_pearson, tmpl=tmpl)
        with mp.Pool(args.threads) as P:
            rv = P.map(apply_pearson_2, d.iterrows())
        r, p = zip(*rv)
        res = pd.DataFrame(index=d.index)
        res['r'] = r
        res['p'] = p
        res['r'] = res['r'].fillna(0)
        res['padj'] = multipletests(res['p'], method='fdr_bh')[1]
        res.sort_values('padj', inplace=True)

    elif args.method == 'euclid':
        res = pd.DataFrame({'r': d.apply(apply_euclid, tmpl=tmpl, axis=1)})
        res.sort_values('padj', inplace=True)

    if binary:
        #add LFC
        lfc_key = res.index.to_series().str.split('__').str.get(0)
        res['junc_lfc'] = list(lfc.loc[lfc_key])
#        res['junc_ttest_t'] = list(norm_tt.loc[lfc_key])
        res['junc_ttest_p'] = list(norm_tp.loc[lfc_key])
        res['junc_ttest_padj'] = list(norm_tpadj.loc[lfc_key])
        res['frac_mean_a'] = meanA
        res['frac_mean_b'] = meanB
        res['frac_diff'] = fracdiff
#        res['frac_ttest_t'] = frac_tt
        res['frac_ttest_p'] = frac_tp
        res['frac_ttest_padj'] = frac_tpadj
    outfile_name = args.template.replace('.template', '') + '.nfj.out'
    lg.info('finished correlation analysis')
    lg.info('writing to: %s', outfile_name)
    res.to_csv(outfile_name, sep="\t", float_format='%.3g')
    lg.info("run time: %.3g seconds", time.time() - start_time)

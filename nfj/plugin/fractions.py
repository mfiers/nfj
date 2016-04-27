#!/usr/bin/env python

import pandas as pd
import logging
from sqlalchemy import create_engine, Index, MetaData 
import leip

lg = logging.getLogger(__name__)

    
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('-n', '--number_to_process', type=int)
@leip.flag('-x', '--exponent', help='convert normalized count values '
           'from log space to count space')
@leip.command
def fraction(app, args):
    """Calculate fractional junction usgae"""

    DEBUG = False
    debug_no = 10000
    
    engine = create_engine(args.datafile)
    lg.info('load from: %s', args.datafile)

    if DEBUG:
        stats = pd.read_sql('select * from junction_stats limit %d' % debug_no,
                            engine, index_col='index')
    else:
        stats = pd.read_sql_table('junction_stats',
                            engine, index_col='index')
    lg.info('loaded %d stats records', len(stats))

    lg.info('load normalized counts')
    if DEBUG:
        d = pd.read_sql('select * from normcounts limit %d' % debug_no,
                        engine, index_col='index')
    else:
        d = pd.read_sql_table('normcounts',
                              engine, index_col='index')


    assert list(stats.index) == list(d.index)

    if args.exponent:
        lg.info("calculating exponent")
        d = 2 ** d

    if d.min().min() < 0:
        lg.error("normalized count table has negative values")
        lg.error("is this log scaled data?? Used voom? If so,")
        lg.error("try running with -x")
        exit()

    assert list(d.index) == list(stats.index)

    def find_forward(r):
        rr = r.split('_')
        rr[2] = 'X'
        return "_".join(rr)

    def find_reverse(r):
        rr = r.split('_')
        rr[1] = 'X'
        return "_".join(rr)

    lg.info('determine forward & reverse aggregation keys')
    stats['forward_group'] = stats.index.to_series().apply(find_forward)
    stats['reverse_group'] = stats.index.to_series().apply(find_reverse)

    assert list(stats.index) == list(d.index)

    # count the number of forward & reverse junctions used for
    # each forward/reverse splice site
    lg.info("count number of forward/reverse junctions")
    fwcount = stats[['forward_group', 'reverse_group']].groupby('forward_group').count()
    fwcount.columns = ['fw_count']
    
    rvcount = stats[['forward_group', 'reverse_group']].groupby('reverse_group').count()
    rvcount.columns = ['rv_count']

    # merge the counts into the splice stat table
    lg.info("merge fw & rv counts into stats table")
    stats = pd.merge(stats, fwcount, how='left', left_on='forward_group', right_index=True)
    stats = pd.merge(stats, rvcount, how='left', left_on='reverse_group', right_index=True)
    lg.info("resulting stats table size: %s", str(stats.shape))

    assert list(stats.index) == list(d.index)
    
    lg.info("sum junctions over all samples / reverse")
    # generate sum for junctions starting or ending on each fw/rv junction
    d_rv = d.copy()
    d_rv['reverse_group'] =  stats['reverse_group']
    
    # one row per unique reverse junction - calculate sum
    d_rv_sum_1 = d_rv.groupby('reverse_group').sum()
    d_rv_sum = pd.merge(stats[['reverse_group']].copy(), d_rv_sum_1,
                        left_on='reverse_group', right_index=True)
    del d_rv_sum['reverse_group']

    lg.info("sum junctions over all samples / forward")
    d_fw = d.copy()
    d_fw['forward_group'] = stats['forward_group']
    
    # one row per unique forward junction - calculate sum
    d_fw_sum_1 = d_fw.groupby('forward_group').sum()
    # expand to one row for each junction
    d_fw_sum = pd.merge(stats[['forward_group']].copy(),
                        d_fw_sum_1, left_on='forward_group', right_index=True)
    del d_fw_sum['forward_group']

    lg.info("calculate fractions")
    d_fw = d.copy()
    d_rv = d.copy()

    d_fw = (d_fw / d_fw_sum).fillna(0)
    d_rv = (d_rv / d_rv_sum).fillna(0)
    
    lg.info("create proper index names, add in gene name")
    d_fw.index = (d_fw.index + '__fw__' + stats['gene'].str.replace(';', '__')).str.rstrip('_')
    stats['forward'] = d_fw.index
  
    d_rv.index = (d_rv.index + '__rv__' + stats['gene'].str.replace(';', '__')).str.rstrip('_')
    stats['reverse'] = d_rv.index
    assert list(stats.index) == list(d.index)
    df = pd.concat([d_fw, d_rv])
    df.index.name = 'index' # for consistency

    lg.info('stats table size: %s', str(stats.shape))
    lg.info('fraction data table size: %s', str(df.shape))

    # remove all non-informative junctions
    dfi = df.copy()

    # for (fw or rv) junction pairs (so not triplets or more) remove one - as one is always
    # 1-x of the other.
    # also - junctions which have no 'pair' are to be removed as well
    fracc_fw = stats[['forward', 'forward_group', 'fw_count']]
    fracc_fw.columns = ['junction', 'group', 'count']
    fracc_rv = stats[['reverse', 'reverse_group', 'rv_count']]
    fracc_rv.columns = ['junction', 'group', 'count']
    
    fracc = pd.concat([fracc_fw, fracc_rv])
    fracc = fracc.sort_values(by='group')

    lg.info("No constitutively used junctions: %d", len(fracc[fracc['count'] == 1]))
    dfi = dfi.loc[fracc[fracc['count'] > 1]['junction']]
    lg.info("After removing them, %d junctions left", len(dfi))

    fracc2 = fracc[fracc['count'] == 2]
    lg.info("Junctions in pairs of two: %d", len(fracc2))
    pair2remove = fracc2.groupby('group').first()['junction']
    lg.info("Removing one half of each pair: %d", len(pair2remove))
    
    dfi = dfi.drop(pair2remove)
    lg.info("After removing, %d junctions left", len(dfi))
    
    # if the total absolute observed variation in the fractional
    # splice site usage, of all samples summed, is less than 1% -
    # remove the junction as uninformative
    # z = (dfi - 1).abs().sum(1)
    # z = ((dfi - dfi.mean()) ** 2).sum(1)
    # print((z > 0.01).value_counts())
    # dfi = dfi[z >= 0.01]
    # lg.info('fraction table after removing uninformative junctions: %s', str(dfi.shape))
    
    if DEBUG: exit()
    if not DEBUG:
        lg.info('start writing fraction data')
        df.to_sql('fraction', engine, if_exists='replace')
        lg.info('start writing informative fraction data')
        dfi.to_sql('fraction_inf', engine, if_exists='replace')
        lg.info('start writing junction stats data')
        stats.to_sql('junction_stats', engine, if_exists='replace')
    else:
        dfi.to_csv('debug_fraction_inf.tsv', sep="\t")
        df.to_csv('debug_fraction.tsv', sep="\t")
        stats.to_csv('debug_fraction_stats.tsv', sep="\t")

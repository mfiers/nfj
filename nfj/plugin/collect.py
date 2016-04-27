
import glob
import logging
import os
import re

import leip
import pandas as pd

from sqlalchemy import create_engine, Index, MetaData 

lg = logging.getLogger(__name__)

GTF = None


@leip.arg('-n', '--nofiles', type=int)
@leip.arg('-N', '--reads_to_process', default=1e10, type=int)
@leip.arg('-r', '--regex', help='regex to identify sample id')
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('indir')
@leip.command
def collect_from_star(app, args):

    engine = create_engine(args.datafile)

    if args.regex:
        re_find_name = re.compile(args.regex)

    def _get_junction_name(row):
        return '%(chrom)s_%(start)s_%(stop)s' % row

    rv = {}

    for i, infile in enumerate(glob.glob('%s/*SJ.out.tab' % args.indir)):
        if args.nofiles and i >= args.nofiles:
            break

        name = os.path.basename(infile).split('SJ')[0]
        lg.info("Processing: %s", name)
        colnames = '''chrom start stop strand intron_motif
                      annotated unique multi max_overhang'''.split()

        d = pd.read_csv(infile, sep="\t", names=colnames)
            
        d.index = d.apply(_get_junction_name, axis=1)
        rv[name] = d['unique']

    lg.info('combine reads from %d input files', len(rv))
    rv = pd.DataFrame(rv).fillna(0)
    lg.info('data shape: %s', str(rv.shape))

    # group by sample
    lg.info('group by sample id')
    rvt = rv.T

    def samplify(f):
        regmatch = re_find_name.search(f)
        sample = regmatch.groups()[0]
        lg.info("filename '%s' -> sample '%s'", f, sample)
        return sample
    
    rvt['sample'] = rvt.index.to_series().apply(samplify)
    
    rvt = rvt.groupby('sample').sum()
    rv = rvt.T.astype(int)
    lg.info("%d junctions observed", len(rv))
    
    # write to disk - raw count table
    lg.info('save to %s', args.datafile)
    rv.to_sql('counts', engine, if_exists='replace')


@leip.arg('-c', '--cutoff', help='fraction has to have at least <cutoff> ' +
          'observations', type=int, default=5)
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.command
def filter_counts(app, args):
    engine = create_engine(args.datafile)

    counts = pd.read_sql('counts', engine, index_col='index')
    lg.info("read %d records", len(counts))

    libsum = counts.sum()
    juncmax = counts.max(1)

    #to junctions per million
    jpm = counts.divide(libsum) * 1e6
    jpm.to_sql('jpm', engine, if_exists='replace')
    
    #remove all which have never more than %d observations
    jpm = jpm[juncmax >= args.cutoff]
    jpm.to_sql('jpm_filter', engine, if_exists='replace')
    lg.info("After filtering: %d jpm records left", len(jpm))
    frac_to_keep = pd.DataFrame(jpm.index.to_series())
    frac_to_keep.to_sql('jpm_keep', engine, if_exists='replace')

    

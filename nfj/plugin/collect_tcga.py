
import glob
import logging
import os

import leip
import pandas as pd

from sqlalchemy import create_engine, Index, MetaData 

lg = logging.getLogger(__name__)

GTF = None


@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('-s', '--sdrf', default='magetab sdrf file')
@leip.arg('indir')
@leip.command
def collect_from_tcga(app, args):
    """
    Collect from tcga precalculated junction count tables
    """

    engine = create_engine(args.datafile)

    rv = {}

    sdrf = pd.read_csv(args.sdrf, sep="\t")
     
    for i, infile in enumerate(glob.glob('%s/*.junction_quantification.txt' % args.indir)):
        sdrf_rec = sdrf[sdrf['Derived Data File'] == os.path.basename(infile)]
        barcode = list(sdrf_rec['Comment [TCGA Barcode]'])[0]
        lg.info("infile %s -> %s", os.path.basename(infile), barcode)

        d = pd.read_csv(infile, sep="\t")
        def fix_junction_name(j):
            a, b = j.split(',')
            chrom, start, strand = a.split(':')
            _, stop, _ = b.split(':')
            start, stop = sorted([int(start), int(stop)])
            return '%s_%d_%d' % (chrom, start, stop)
        
        d['junction'] = d['junction'].apply(fix_junction_name)
        d = d.groupby('junction').sum()['raw_counts']
        rv[barcode] = d
            
    lg.info('combine reads from %d input files', len(rv))
    rv = pd.DataFrame(rv).fillna(0)
    rv.index.name = 'index'

    lg.info('data shape: %s', str(rv.shape))

    lg.info("remove junctions with total count < 5")
    rv = rv[rv.sum(1) > 4]

    lg.info('data shape: %s', str(rv.shape))

    lg.info("%d junctions observed", len(rv))
    
    # write to disk - raw count table
    lg.info('save to %s', args.datafile)
    rv.to_sql('counts', engine, if_exists='replace')


# @leip.arg('-c', '--cutoff', help='fraction has to have at least <cutoff> ' +
#           'observations', type=int, default=5)
# @leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
# @leip.command
# def filter_counts(app, args):
#     engine = create_engine(args.datafile)

#     counts = pd.read_sql('counts', engine, index_col='index')
#     lg.info("read %d records", len(counts))

#     libsum = counts.sum()
#     juncmax = counts.max(1)

#     #to junctions per million
#     jpm = counts.divide(libsum) * 1e6
#     jpm.to_sql('jpm', engine, if_exists='replace')
    
#     #remove all which have never more than %d observations
#     jpm = jpm[juncmax >= args.cutoff]
#     jpm.to_sql('jpm_filter', engine, if_exists='replace')
#     lg.info("After filtering: %d jpm records left", len(jpm))
#     frac_to_keep = pd.DataFrame(jpm.index.to_series())
#     frac_to_keep.to_sql('jpm_keep', engine, if_exists='replace')

    


import logging
import leip
import pandas as pd

from nfj import util

lg = logging.getLogger(__name__)

@leip.arg('--db', default='nfj')
@leip.command
def stats(app, args):
    """Calculate stats on the raw & normalized counts table"""

    counts= util.load(args.db, 'counts')
    no_samples = counts.shape[1]

    lg.info("read %d records", len(counts))
    cstats = pd.DataFrame(index=counts.index)

    coords = cstats.index.str.rsplit('_',2)

    lg.info('parse junction name')
    cstats['chrom'] = coords.str.get(0)
    cstats['start'] = coords.str.get(1).astype(int)
    cstats['stop'] = coords.str.get(2).astype(int)
    
    lg.info('raw counts: calculate mean')
    cstats['counts_mean'] = counts.mean(1)

    lg.info('raw counts: calculate min')
    cstats['counts_min'] = counts.min(1)
    
    lg.info('raw counts: calculate max')
    cstats['counts_max'] = counts.max(1)

    lg.info('raw counts: standard deviation')
    cstats['counts_std'] = counts.std(1)

    lg.info('raw counts: coefficient of variation')
    cstats['counts_cv'] = cstats['counts_std']  / cstats['counts_mean']

#    for h in [5, 10]:
#        lg.info('raw counts: calculate high_%02d' % h)
#        if no_samples < h-1:
#            break
#        cstats['counts_high_%02d' % h] = counts.apply(lambda x: x.sort_values()[-h], axis=1)

    lg.info("read normalized count table")
    normcounts= util.load(args.db, 'normcounts')

    lg.info("normalized counts: read %d records", len(normcounts))
    
    lg.info('normalized counts: calculate mean')
    cstats['norm_mean'] = normcounts.mean(1)
    
    lg.info('normalized counts: calculate min')
    cstats['norm_min'] = normcounts.min(1)
    
    lg.info('normalized counts: calculate min')
    cstats['norm_max'] = normcounts.max(1)
    
    lg.info('normalized counts: standard deviation')
    cstats['norm_std'] = normcounts.std(1)

    lg.info('normalized counts: coefficient of variation')
    cstats['norm_cv'] = cstats['norm_std']  / cstats['norm_mean']

    util.save(args.db, junction_stats=cstats)

    
RSCRIPT="""
library(edgeR)
library(limma)

d = read.table('_tmp_voom_input.tsv', sep="\t", header=1, row.names=1)
dge=DGEList(counts=d)
dge = calcNormFactors(dge)
v = voom(dge, design=NULL, plot=FALSE)
write.table(v$E, '_tmp_voom_output.tsv', sep="\t")
"""
    
@leip.arg('--db', default='nfj')
@leip.command
def voom(app, args):
    "Normalize using Limma's voom"
    from sh import Rscript

    counts= util.load(args.db, 'counts')
    lg.info("read %d records", len(counts))
    
    lg.info("write raw counts to tmp file")
    counts.to_csv('_tmp_voom_input.tsv', sep="\t")

    lg.info("start R, normalize using voom")
    Rscript('-', _in=RSCRIPT)

    lg.info("read R output")
    v = pd.read_csv('_tmp_voom_output.tsv', sep="\t")

    #back to count space (from log space)
    v = 2 ** v
    
    lg.info("Write to the database")
    util.save(args.db, normcounts=v)

    
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.command
def jpm(app, args):
    "Normalize to junctions per million (like RPM)"
    lg.critical("not implemented")
    exit(-1)
    engine = create_engine(args.datafile)
    counts = pd.read_sql('counts', engine, index_col='index')
    lg.info("read %d records", len(counts))
    
    #to junctions per million
    jpm = counts.divide(libsum) * 1e6
    lg.info("Write to the database")
    jpm.to_sql('normcounts', engine, if_exists='replace')

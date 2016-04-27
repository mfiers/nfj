
import logging
import warnings

import leip
import pandas as pd

warnings.simplefilter(action="ignore", category=FutureWarning)

lg = logging.getLogger(__name__)


@leip.arg('-f', '--no_flank', default=10, type=int)
@leip.arg('-d', '--datafile', default='nfj.h5')
@leip.arg('fjunction')
@leip.command
def fplot(app, args):

    # late import - makes the system a little quicker
    import matplotlib.pyplot as plt
    import seaborn as sns

    jnc = args.fjunction
    lg.info('junction: %s', jnc)

    chrom, start, stop, _ = jnc.split('_', 3)
    start, stop = int(start), int(stop)
    
    lg.info("chromosome: %s" % chrom)
    lg.info("start: %d" % start)
    lg.info("stop: %d" % stop)

    query1 = "chr = %r & start > %d & stop < %d" % (
        chrom, start-1e6, stop+1e6)
    lg.info('get flanking: %s', query1)
        
    with pd.HDFStore(args.datafile, mode='r') as store:
       d = store.select('sstats', query1)

    assert len(d) > 0
    d = d.sort_values(by=['start', 'stop'])
    sjnc = "_".join(jnc.split('_')[:4])

    assert sjnc in d.index
    spos = list(d.index).index(sjnc)
    ssta = max(0, spos - args.no_flank)
    ssto = min(len(d.index)-1, spos+args.no_flank)

    d = d.iloc[ssta:ssto]
    print(d.iloc[:, :])
    print(spos)
    exit()
#    flanking = 
    
    query = 'index = %r' % args.fjunction
    lg.info('query: %s' % query)

    #find nearby indici
#        exit()
    
    with pd.HDFStore(args.datafile, mode='r') as store:
        d = store.select('fraction', where=query)

    lg.info("result shape: %s", d.shape)
    assert d.shape[0] == 1
    
    d = d.T
    d.columns = ['frac']
    d['group'] = d.index.str[:8]

    sns.boxplot(x='group', y='frac', data=d)
    sns.stripplot(x="group", y="frac", data=d, size=5,
                  jitter=True, edgecolor='k')
     
    plt.title('junction fraction usage plot for: %s' % args.fjunction)
    ymin, ymax = plt.ylim()
    if ymin < 0.05:
        ymin = -0.05
    if ymax > 0.95:
        ymax = 1.05
        
    plt.ylim(ymin, ymax)
    
    plotname = 'fplot_%s.png' % args.fjunction
    plt.savefig(plotname)
    
    

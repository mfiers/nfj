
import logging

from path import Path
import leip

lg = logging.getLogger(__name__)

    
RSCRIPT="""
library(edgeR)
library(limma)

d = read.table('featurecounts.grouped.tsv', sep="\t", header=1, row.names=1)
dge=DGEList(counts=d)
dge = calcNormFactors(dge)
v = voom(dge, design=NULL, plot=FALSE)
write.table(v$E, 'featurecounts.norm.tsv', sep="\t")
"""

@leip.arg('--db', default='./nfj')
@leip.arg('--gtfdb', default='gtf')
@leip.arg('-r', '--regex')
@leip.arg('indir')
@leip.command
def prep_gene_expr_star(app, args):
    """
    Perform differential gene expression
    """

    import re
    import pandas as pd
    
    from sh import Rscript
    from sh import featureCounts
    from nfj import util

    gtfdir = Path(args.gtfdb)

    fc_out = Path('featurecounts.out.tsv')
    if not fc_out.exists():
        lg.info("Running featurecounts")
        featureCounts(
            '-a', gtfdir / 'full.gtf',
            '-g', 'gene_name', '-Q', 10,
            '-T', 20, '-o', 'featurecounts.out.tsv',
            *Path(args.indir).glob('*Aligned.sortedByCoord.out.bam'))


    if not Path("featurecounts.grouped.tsv").exists():
        lg.info("read raw gene counts")

        d = pd.read_csv('featurecounts.out.tsv', sep="\t", comment='#',
                        index_col=0)
        d.index.name = 'gene'
        dmeta = d.iloc[:,:5]
        dcount = d.iloc[:,6:]
        _t = dcount.T
        rex = re.compile(args.regex)
        _t['sample'] = [rex.search(Path(x).basename()).groups()[0]
                        for x in _t.index]
        dcount = _t.groupby('sample').sum().T
        dcount.to_csv('featurecounts.grouped.tsv', sep="\t")

        util.save(args.db,
                  genecounts_meta = dmeta,
                  genecounts_raw = dcount)


    if not Path("featurecounts.norm.tsv").exists():
        lg.info("start R, normalize using voom")
        Rscript('-', _in=RSCRIPT)
    
    lg.info("read R output")
    v = pd.read_csv("featurecounts.norm.tsv", sep="\t")
    print(v.iloc[:5,:5])

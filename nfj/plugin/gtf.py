
from collections import defaultdict
import logging
lg = logging.getLogger(__name__)


import pandas as pd
import numpy as np
import leip
from sqlalchemy import create_engine, Index, MetaData 

@leip.flag('--debug')
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('gtf')
@leip.command
def load_gtf(app, args):

    engine = create_engine(args.datafile)
    lg.info('loading: %s' % args.gtf)
    
    gtfcols = 'chr source type start stop score strand frame attrs'.split()
    gtfargs = dict(sep= "\t", names= gtfcols, comment='#')
    gtfargs['dtype'] = {'start': np.uint32, 'stop': np.uint32, 'chr': str}
    if args.debug:
        gtf = pd.read_csv(args.gtf, nrows=1000, **gtfargs)
    else:
        gtf = pd.read_csv(args.gtf, **gtfargs)

    lg.info('loaded %d records' % len(gtf))
    genes = gtf[gtf['type'] == 'gene']
    transcripts = gtf[gtf['type'] == 'transcript']
    exons = gtf[gtf['type'] == 'exon']
    
    def _find_thing(x, what='gene_name'):
        x = x.split(what)[1].strip().split('"')[1]
        return x
    
    genes['gene_name'] = genes['attrs'].apply(_find_thing, what='gene_name')
    genes['gene_id'] = genes['attrs'].apply(_find_thing, what='gene_id')
    genes['biotype'] = genes['attrs'].apply(_find_thing, what='gene_biotype')
    del genes['attrs']
    
    exons['gene_id'] = exons['attrs'].apply(_find_thing, what='gene_id')
    exons['transcript_name'] = exons['attrs'].apply(_find_thing, what='transcript_name')
    exons['transcript_biotype'] = exons['attrs'].apply(_find_thing, what='transcript_biotype')
    exons['transcript_id'] = exons['attrs'].apply(_find_thing, what='transcript_id')
    exons['gene_name'] = exons['attrs'].apply(_find_thing, what='gene_name')
    exons['gene_biotype'] = exons['attrs'].apply(_find_thing, what='gene_biotype')
    exons['gene_id'] = exons['attrs'].apply(_find_thing, what='gene_id')
    del exons['attrs']
    
    # generate a per gene collection of exons, junctions & biotypes
    junction_effect_table = []
    
    for _, generec in genes.iterrows():
        gene_id = generec['gene_id']
        gene_name = generec['gene_name']
        xtable = {}
        jtable = {}
        btypes = defaultdict(lambda: set())
        xt = exons[exons['gene_id'] == gene_id]
        for transid in xt['transcript_id'].unique():
            xtt = xt[xt['transcript_id'] == transid]
            if len(xtt) == 1:
                #one exon transcript - no junctions - not interseting
                continue
            biotype = list(xtt['transcript_biotype'].unique())                
            assert len(biotype) == 1
            biotype = biotype[0]
            xtt = xtt.sort_values(by=['chr', 'start', 'stop'])
            exon_coords = [(x['chr'], x['start'], x['stop'])
                           for (nx, x) in xtt.iterrows()]
            exon_sets = list(zip(exon_coords[:-1], exon_coords[1:]))
            
            assert min([x[1][1] - x[0][2] for x in exon_sets]) >= 0
            xtable[transid] = exon_coords
            jtable[transid] = exon_sets
            btypes[biotype].update(set(exon_coords))

        #invert - exon->biotype table
        exon2biotype = defaultdict(set)        
        for b in btypes:
            for x in btypes[b]:
                exon2biotype[x].add(b)

        # now associate unique biotypes & calculate exon naems
        for transid in jtable:
            for jstart, jstop in jtable[transid]:
                bt_start = exon2biotype[jstart]
                bt_stop = exon2biotype[jstop]
                if len(bt_start) == 1 and len(bt_stop) == 1:
                    # I assume this is the case - but if I'm mistaken...
                    assert bt_start == bt_stop
                    
                if len(bt_start) == 1:
                    jat = list(bt_start)[0]
                elif len(bt_stop) == 1:
                    jat = list(bt_stop)[0]
                else:
                    continue
                jname = '%s_%s_%s' % (jstart[0], jstart[2]+1, jstop[1]-1)
                junction_effect_table.append([gene_id, gene_name, jname, jat])
                
    junction_effect_table = pd.DataFrame.from_records(junction_effect_table,
                                         columns='gene_id gene_name junction_name effect'.split())

    junction_effect_table.to_csv('test.tsv', sep="\t")
    if not args.debug:
        junction_effect_table.to_sql('gtf_junction_effect', engine, if_exists='replace')
        exon.to_sql('gtf_exon', engine, if_exists='replace')
        genes.to_sql('gtf_genes', engine, if_exists='replace')
        
        #add indici
        meta = MetaData()
        meta.reflect(engine)
        table = meta.tables['gtf_genes']
        i = Index('genes.pos', table.c.chr, table.c.start, table.c.stop)
        i.create(engine)

    
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.command
def find_gene(app, args):

    lg.info("start junction annotation")
    engine = create_engine(args.datafile)

    stats = pd.read_sql('SELECT * FROM junction_stats',
                        engine, index_col='index')

    lg.info("loaded %d records ", len(stats))

    gtf = pd.read_sql('gtf', engine)
    
    lg.info("loaded %d records from the gtf table", len(gtf))

    from multiprocessing.dummy import Pool
    def _find_gene(namerow):
        name, row = namerow
        chrom, start, stop = row['chrom'], row['start'], row['stop']

        start, stop = list(sorted([start, stop]))
        grec = gtf[(chrom == gtf['chr']) &
                   (start < gtf['stop']) &
                   (stop > gtf['start'])]
        if len(grec) > 0:
            return ";".join(grec['name'])
        else:
            return ""

    with Pool(20) as p:
        rv = p.map(_find_gene, stats.iterrows())

    stats['gene'] = rv
    lg.info("update junction_stats table")
    stats.to_sql('junction_stats', engine, if_exists='replace')

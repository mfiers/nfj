
from collections import defaultdict
import logging
lg = logging.getLogger(__name__)

import leip
from path import Path
from nfj import util

@leip.flag('--debug')
@leip.arg('--gtfdb', default='gtf')
@leip.arg('gtf')
@leip.command
def load_gtf(app, args):

    import pandas as pd
    import numpy as np

    lg.info('loading: %s' % args.gtf)
    
    gtfcols = 'chr source type start stop score strand frame attrs'.split()
    gtfargs = dict(sep= "\t", names= gtfcols, comment='#')
    gtfargs['dtype'] = {'start': np.uint32, 'stop': np.uint32, 'chr': str}
    
    if args.debug:
        gtf = pd.read_csv(args.gtf, nrows=5000, **gtfargs)
    else:
        gtf = pd.read_csv(args.gtf, **gtfargs)

    lg.info('loaded %d records' % len(gtf))
    
    def _find_thing(x, what='gene_name'):
        x = x.split(what)[1].strip().split('"')[1]
        return x

    # 1_10025064_10027102     (('ENSMUSG00000025917', 'Cops5', 'nonsense_mediated_decay'),
    #            ('ENSMUSG00000025917', 'Cops5', 'protein_coding'))

    lg.info('create gene table')
    genes = gtf[gtf['type'] == 'gene'].copy()
    genes['gene_name'] = genes['attrs'].apply(_find_thing, what='gene_name')
    genes['gene_id'] = genes['attrs'].apply(_find_thing, what='gene_id')
    genes['biotype'] = genes['attrs'].apply(_find_thing, what='gene_biotype')
    genes.sort_values(by='gene_id', inplace=True)
    lg.info('found %d genes', len(genes))
    del genes['attrs']

    lg.info('create exon table')
    exons = gtf[gtf['type'] == 'exon'].copy()
    exons['gene_id'] = exons['attrs'].apply(_find_thing, what='gene_id')
    exons['transcript_name'] = exons['attrs'].apply(_find_thing, what='transcript_name')
    exons['transcript_biotype'] = exons['attrs'].apply(_find_thing, what='transcript_biotype')
    exons['transcript_id'] = exons['attrs'].apply(_find_thing, what='transcript_id')
    exons['gene_name'] = exons['attrs'].apply(_find_thing, what='gene_name')
    exons['gene_biotype'] = exons['attrs'].apply(_find_thing, what='gene_biotype')
    exons['gene_id'] = exons['attrs'].apply(_find_thing, what='gene_id')
    lg.info('found %d exons', len(exons))
    del exons['attrs']
    
    # generate a per gene collection of exons, junctions & biotypes
    junction_effect_table = set()
    lg.info('create a exon biotype table')
    no_genes = len(genes)
    for i, (_, generec) in enumerate(genes.iterrows()):
        gene_id = generec['gene_id']
        gene_name = generec['gene_name']
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
            jtable[transid] = exon_sets
            btypes[biotype].update(set(exon_coords))

        #invert - exon->biotype table
        exon2biotype = defaultdict(set)        
        for b in btypes:
            for x in btypes[b]:
                exon2biotype[x].add(b)

        # now go from exon2biotype to junction2biotype
        jstart2biotype = defaultdict(set)
        jstop2biotype = defaultdict(set)


        for transid in jtable:
            for jstart, jstop in jtable[transid]:
                
                bt_start = exon2biotype[jstart]
                bt_stop = exon2biotype[jstop]
               
                jstart2biotype[jstart[2]+1].update(bt_start)
                jstop2biotype[jstop[1]-1].update(bt_stop)
             
        for transid in jtable:
            for jstart, jstop in jtable[transid]:

                jname = '%s_%s_%s' % (jstart[0], jstart[2]+1, jstop[1]-1)
                bt_start = jstart2biotype[jstart[2]+1]
                bt_stop = jstop2biotype[jstop[1]-1]
                
                if len(bt_start) == 1 and len(bt_stop) == 1:
                    # I assume this is the case - but if I'm mistaken...
                    assert bt_start == bt_stop
                    
                if len(bt_start) == 1:
                    jat = list(bt_start)[0]
                elif len(bt_stop) == 1:
                    jat = list(bt_stop)[0]
                else:
                    jat = 'not_informative'

                if jname == '1_10047527_10058719':
                    print(jname, gene_id, gene_name, jat)
                    print(jstart, bt_start)
                    print(jstop, bt_stop)
                junction_effect_table.add((gene_id, gene_name, jname, jat))

    junction_effect_table = pd.DataFrame.from_records(list(junction_effect_table),
                                         columns='gene_id gene_name junction_name effect'.split())
    lg.info("found %d junction effects", len(junction_effect_table))

    util.save(args.gtfdb,
              junction_effect=junction_effect_table,
              exon=exons,
              gene=genes)


def _find_effect(group):
    name = group[0]
    group = group[1].drop_duplicates()
    return name, tuple(set(group.apply(lambda r: (r['gene_id'], r['gene_name'], r['effect']), axis=1)))


def _find_gene(namerow, genetable, flank=0, tag='_overlap_'):
    import numpy as np
    name, row = namerow
    chrom, start, stop = row['chrom'], row['start']-flank, row['stop']+flank

    start, stop = list(sorted([start, stop]))
    grec = genetable[(chrom == genetable['chr']) &
               (start < genetable['stop']) &
               (stop > genetable['start'])]
    if len(grec) == 0:
        return np.nan
    rv = grec.apply(lambda x: (x['gene_id'], x['gene_name'], tag), axis=1)
    return tuple(set(rv))



@leip.flag('--ucsc', help='fix chromosome & junction names to use ensembl annotation with a ucsc mapping')
@leip.flag('--debug')
@leip.arg('--threads', '-j', type=int, default=20)
@leip.arg('--gtfdb', default='gtf')
@leip.arg('--db', default='nfj')
@leip.command
def annotate(app, args):

    from functools import partial
    from multiprocessing import Pool
    
    import pandas as pd
    import numpy as np

    lg.info("start junction annotation")

    stats= util.load(args.db, 'junction_stats')

    lg.info("loaded %d observed junction records", len(stats))

    gtf_genes = util.load(args.gtfdb, 'gene')
    gtf_genes = gtf_genes.sort_values(by='chr start stop'.split())
    lg.info("loaded %d gene records", len(gtf_genes))

    junction = util.load(args.gtfdb, 'junction_effect')
    lg.info("loaded %d junction records", len(junction))

    if args.ucsc:
        lg.info("prepending chr to chromosome names")
        gtf_genes['chr'] = 'chr' + gtf_genes['chr']
        junction['junction_name'] = 'chr' + junction['junction_name']

        # subtract one from all annotation tables :(
        gtf_genes['start'] = gtf_genes['start'] - 1
        gtf_genes['stop'] = gtf_genes['stop'] - 1
        def _jnf(n):
            a,b,c = n.rsplit('_',2)
            return '%s_%d_%d' % (a, int(b)-1, int(c)+1)
        junction['junction_name'] = junction['junction_name'].apply(_jnf)
    junction = junction.sort_values(by='junction_name')
    no_unqiue_junc = len(set(junction['junction_name']))
    lg.info(" -- containing %d unique junctions", no_unqiue_junc)

    annotated = pd.DataFrame(index=stats.index)

    overlap = set(stats.index) & set(junction['junction_name'])
    lg.info("found %d (%.2f%%) observed junctions in annotation",
            len(overlap), 100 * (len(overlap) / no_unqiue_junc))

    junction = junction.sort_values(by='junction_name')

    #group predicted junction effects on junction
    junction = junction[junction['junction_name'].isin(overlap)]

    lg.info("after pruning unobserved junctions: %d left", len(junction))
    if args.debug:
        jgroup = junction.head(1000).groupby('junction_name')
    else:
        jgroup = junction.groupby('junction_name')
   
    lg.info("start link annotated to observed junctions")
    with Pool(args.threads) as P:
        effect = P.map(_find_effect, jgroup)
        lg.info("found %d junction effects" % len(effect))

    effect_index, effect = zip(*effect)
    effect = pd.Series(effect, index=effect_index)
    annotated['effect'] = effect

    no_ann = len(annotated.dropna())
    lg.info("total of %d junctions annotated, (%d to go)",
            no_ann, len(annotated) - no_ann)
    
    not_annotated = annotated[pd.isnull(annotated['effect'])]
    lg.info("finding gene overlap for %d unannotated junctions", len(not_annotated))
    if args.debug:
        not_annotated = not_annotated.head(100)
        lg.info("debug mode - only doing first 100")
        
    with Pool(args.threads) as P:
        find_gene_2 = partial(_find_gene, genetable=gtf_genes, tag='_overlap_')
        gene_2 = P.map(find_gene_2, stats.loc[not_annotated.index].iterrows())

    gene_2 = pd.Series(gene_2, index=not_annotated.index)
    lg.info("discovered %d new overlaps", len(gene_2.dropna()))
    annotated.loc[gene_2.index, 'effect'] = gene_2
       
    #now again, but within 5k of a gene
    not_annotated = annotated[pd.isnull(annotated['effect'])]
    lg.info("finding gene 5k proximity for %d unannotated junctions", len(not_annotated))
    if args.debug:
        not_annotated = not_annotated.head(100)
        lg.info("debug mode - only doing first 100")

    with Pool(args.threads) as P:
        find_gene_2 = partial(_find_gene, genetable=gtf_genes, flank=5000, tag='_flank_5k_')
        gene_2 = P.map(find_gene_2, stats.loc[not_annotated.index].iterrows())
        
    gene_2 = pd.Series(gene_2, index=not_annotated.index)
    lg.info("discovered %d new flanking", len(gene_2.dropna()))
    annotated.loc[gene_2.index, 'effect'] = gene_2

    lg.info("update junction_stats table")
    util.save(args.db, junction_annotation=annotated)
#    stats.to_sql('junction_stats', engine, if_exists='replace')


import logging
import warnings

import leip

from nfj import util
warnings.simplefilter(action="ignore", category=FutureWarning)

lg = logging.getLogger(__name__)

def item_overlap(item1, item2, spacer=0):
    rv = (item1[1] + spacer) >= item2[0] and (item1[0] - spacer) <= item2[1]
    return rv
def item_overlap(item1, item2, spacer=0):
    rv = (item1[1] + spacer) >= item2[0] and (item1[0] - spacer) <= item2[1]
    return rv

def rowify(rowdb, new_item, spacer=2):
    for row_no, row in enumerate(rowdb):
        for row_item in row:
            if item_overlap(row_item, new_item, spacer): break
        else: # no overlap - add
            row.append(new_item)
            return row_no
    else:
        rowdb.append([new_item])
        return len(rowdb)-1

@leip.arg('--db', default='nfj')
@leip.arg('--gtfdb', default='gtf')
@leip.flag('--ucsc')
@leip.arg('--name', nargs=2)
@leip.arg('gene')
@leip.arg('template')
@leip.command
def eggplot(app, args):

    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    import matplotlib.gridspec as gridspec
    from matplotlib.path import Path

    import seaborn as sns

    names = args.name

    gene = args.gene
    rowdb = []

    template = pd.read_csv(args.template, header=None, sep="\t",
                           index_col=0).iloc[:,0].astype(int)

    stats = util.load(args.db, 'stats2')
    lg.info("loaded stats: %s", stats.shape)
    ncounts = util.load(args.db, 'normcounts')
    lg.info("loaded normalized counts: %s", ncounts.shape)
    fraction = util.load(args.db, 'fraction')
    lg.info("loaded fractions: %s", fraction.shape)

    genedb = util.load(args.gtfdb, 'gene')
    lg.info("loaded genes: %s", genedb.shape)
    sgene = genedb[genedb['gene_name'] == gene].iloc[0]
    lg.info('%s gene found: chrom %s: %d-%d', sgene['gene_name'],
            sgene['chr'], sgene['start'], sgene['stop'])
    lg.info('ensembl id: %s', sgene['gene_id'])

    if args.ucsc:
        sgene['chr'] = 'chr' + sgene['chr']

    flank = 50000
    sstat = stats[ (stats['chrom'] == sgene['chr']) &
                   (stats['start'] < sgene['stop'] + flank) &
                   (stats['stop'] > sgene['start'] - flank) ]


    lg.info("stats records found: %d", len(sstat))

    snorm = ncounts.loc[sstat.index]
    sstat = sstat.loc[snorm.index]
    fracnames = pd.concat([sstat['forward'], sstat['reverse']])
    sfrac = fraction[fraction.index.to_series().isin(fracnames)]
    lg.info("no select records in fraction table: %d", len(sfrac))

    plt.figure(figsize=(14,5))
    gs = gridspec.GridSpec(2, 1, hspace=0.05, height_ratios=[2,2])
    ax_nd = plt.subplot(gs[0,:])
    ax_fp = plt.subplot(gs[1,:])

    allX = set()

    gvc = template.value_counts()
    groups = list(reversed(sorted(gvc.index)))
    assert len(groups) == 2
    template = template.sort_index()
    template = template.sort_values(ascending=False)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        palette = sns.color_palette(n_colors=7)

    if True: # calculate normalized * fraction mean/sd
        snmean = pd.DataFrame(index=snorm.index)
        snorm_select = snorm[template.index]
        for g in groups:
            sn_sg = snorm_select.loc[:,template == g]
            snmean['mean_%s' % g] = sn_sg.mean(1)
            snmean['sd_%s' % g] = sn_sg.std(1)
        fmean = pd.DataFrame(index=sfrac.index)
        sfrac_select = sfrac[template.index]
        for g in groups:
            sn_sg = sfrac_select.loc[:,template == g]
            fmean['mean_%s' % g] = sn_sg.mean(1)
            fmean['sd_%s' % g] = sn_sg.std(1)
    if True: # Norm Junction Count Plot
        curveCode = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        curveCode2 = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                      Path.CURVE4, Path.CURVE4, Path.CURVE4]
        allmean = snmean[['mean_%s' % x for x in groups]]

        mxj = allmean.max().max(); mnj =allmean.min().min()
        g1, g2 = groups
        mean1 = snmean['mean_%s' % g1]
        mean2 = snmean['mean_%s' % g2]
        sd1 = snmean['sd_%s' % g1]
        sd2 = snmean['sd_%s' % g2]
        label_seen = set()

        for i2, (n1, j1) in enumerate(mean1.iteritems()):
            j2 = mean2[n1]
            srec = sstat.loc[n1]
            start, stop = sstat.loc[n1, 'start'], sstat.loc[n1, 'stop']
            allX.add(start); allX.add(stop)
            ptchpar2 = {}
            if j1 > j2:
                fillcolor = palette[0]
                if not 'a' in label_seen:
                    ptchpar2 = dict(label='more in %s' % names[g1])
                    label_seen.add('a')
            else:
                fillcolor = palette[1]
                if not 'b' in label_seen:
                    ptchpar2 = dict(label='more in %s' % names[g2])
                    label_seen.add('b')
            pthd = Path([(start, mnj), (start, j1), (stop, j1), (stop,mnj),
                         (stop, j2), (start, j2), (start,mnj)], curveCode2)
            ax_nd.add_patch(patches.PathPatch(pthd, lw=1, ec=fillcolor, fc=fillcolor,
                                              alpha=1, **ptchpar2))
            mxj = max(mxj, j1, j2)
        ax_nd.set_ylabel("Norm.Jnc.Cnt")
        ax_nd.legend(loc='upper center', ncol=2)
        ax_nd.set_xlim(min(allX), max(allX))
        ax_nd.set_ylim(mnj, mxj*1.3)

    if True: # Fraction plot
        assert len(groups) == 2
        palette = sns.color_palette(n_colors=len(groups))
        curveCode = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO,
                     Path.CURVE3, Path.CURVE3]
        curveFH = [Path.MOVETO, Path.CURVE3, Path.CURVE3]

        means = pd.DataFrame(index=sfrac.index)
        stds =  pd.DataFrame(index=sfrac.index)
        maxdiff = 0
        allY = []

        for j, (jname, jstat) in enumerate(sstat.iterrows()):
            ffrac = sfrac.loc[jstat['forward']]
            rfrac = sfrac.loc[jstat['reverse']]
            start, stop = jstat[['start', 'stop']]
            halfway = (0.5 * (stop - start)) + start
            fmeans = pd.Series(index=groups)
            fstds = pd.Series(index=groups)
            rmeans = pd.Series(index=groups)
            rstds = pd.Series(index=groups)
            for i, group in enumerate(groups):
                samples = template[template == group].index
                fmeans[group] = ffrac[samples].mean()
                fstds[group] = ffrac[samples].std()
                rmeans[group] = rfrac[samples].mean()
                rstds[group] = rfrac[samples].std()

            ddmeans = pd.concat([fmeans, rmeans], axis=1)
            ddmeans.columns = ['f', 'r']
            dmeans = ddmeans.loc[groups[1]] - ddmeans.loc[groups[0]]
            maxdiff = max(dmeans.abs().sum(), maxdiff)
            df = dmeans['f']; dr = dmeans['r']

            if not isinstance(jstat['effect'], (list, tuple)):
                effext = '?'
            elif len(jstat['effect']) > 1:
                efftext = str(jstat['effect'])
            elif jstat['effect'][0][1] == gene:
                efftext = jstat['effect'][0][2]
            else:
                efftext = jstat['effect'][0][1] + "_" + jstat['effect'][0][2]

            print(df, dr, efftext)
            if abs(df) > 0.02:
                pth = Path([(start, 0), (start, df), (halfway, df)], curveFH)
                allY.append(df)
                fillcolor = palette[1] if df > 0 else palette[0]
                ax_fp.add_patch(patches.PathPatch(pth, lw=4, ec=fillcolor, fc='none', alpha=1))
                pth = Path([(stop, 0), (stop, df), (halfway, df)], curveFH)
                ax_fp.add_patch(patches.PathPatch(pth, lw=3, ec=fillcolor,  ls=':', fc='none', alpha=0.8))
                va = 'bottom' if df > 0 else 'top'
                ax_fp.text(halfway, df, efftext, ha='center', va=va)
            if abs(dr) > 0.02:
                pth = Path([(stop, 0), (stop, dr), (halfway, dr)], curveFH)
                allY.append(dr)
                fillcolor = palette[1] if dr > 0 else palette[0]
                ax_fp.add_patch(patches.PathPatch(pth, lw=4, ec=fillcolor, fc='none', alpha=1))
                pth = Path([(start, 0), (start, dr), (halfway, dr)], curveFH)
                ax_fp.add_patch(patches.PathPatch(pth, lw=3, ec=fillcolor, ls=':', fc='none', alpha=0.8))
                va = 'bottom' if dr > 0 else 'top'
                ax_fp.text(halfway, dr, efftext, ha='center', va=va)

        if allY:
            ax_fp.set_ylim(min(allY)*1.1, max(allY)*1.1)
        ax_fp.set_xlim(min(allX), max(allX))
        ax_fp.set_ylabel("dFJ")
        ax_fp.axhline(0, color='grey')
    plt.savefig(gene + '.png', dpi=200)

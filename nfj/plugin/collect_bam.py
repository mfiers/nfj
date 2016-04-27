
from collections import Counter
import glob
import logging
import os
import re

import leip

import pandas as pd
import pysam
from path import Path

from sqlalchemy import create_engine, Index, MetaData 

lg = logging.getLogger(__name__)

GTF = None

@leip.arg('-n', '--nofiles', type=int)
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('-N', '--reads_to_process', default=1e10, type=int)
@leip.arg('-q', '--mapping_quality', default=10, type=int)
@leip.arg('-r', '--regex', help='regex to identify sample id')
@leip.arg('indir')
@leip.command
def collect_from_bam(app, args):

    if args.regex:
        re_find_name = re.compile(args.regex)

        
    engine = create_engine(args.datafile)
    def _get_junction_name(row):
        return '%(chrom)s_%(start)s_%(stop)s_%(strand)s' % row

    junctions = {}

    for i, infile in enumerate(Path(args.indir).expanduser().glob('*.bam')):
        
        lg.info("processing bam file %s", infile)

        name = infile.basename()

        if args.nofiles and i >= args.nofiles:
            break

        bamfile = pysam.AlignmentFile(infile, "rb")
        jtable = Counter()
        for j, read in enumerate(bamfile):
            if j > args.reads_to_process:
                break
            if read.mapping_quality < args.mapping_quality:
                continue
            cigevents = [x[0] for x in read.cigartuples]
            no_junctions = cigevents.count(3)
            if no_junctions < 1:
                continue
            pairs = dict(read.get_aligned_pairs(matches_only=True))

            pos_in_read = 0
            junctions_seen = 0
            primed = False
            for event, rpos in read.cigartuples:
                if event in [0, 1, 4]:
                    pos_in_read += rpos
                    continue
                elif event in [2]:
                    # deletion
                    primed = True
                    continue
                elif event == 3:
                    junctions_seen += 1
                    jstart = pairs[pos_in_read-1]+2
                    jstop =  pairs[pos_in_read]
                    jname = '%s_%d_%d' % (read.reference_name, jstart, jstop)
                    if primed in jtable:
                        lg.critical("primed junction: " + jname)
                        exit()
#                    if not jname in jtable:
#                        lg.info("Found junction %s", jname)
                    jtable[jname] += 1
                else:
                    lg.critical("unkonwn event: %d %d", event, rpos)
                    print(j, pos_in_read, read.reference_name)
                    print(pairs)
                    print(read)
                    print(read.cigarstring)
                    print(read.cigartuples)
                    exit()
                if junctions_seen >= no_junctions:
                    break

        col = pd.Series(jtable)
        lg.info("%d reads processed", j+1)
        lg.info("%d unique junctions observed", col.shape[0]+1)
        lg.info("%d reads covering all junctions", col.sum())
        junctions[name] = col

    lg.info('combine into one dataframe from %d input files', len(junctions))
    rv = pd.DataFrame(junctions).fillna(0)
    lg.info('data shape: %s', str(rv.shape))

    # group by sample
    lg.info('group by sample id')
    rvt = rv.T

    def samplify(f):
        regmatch = re_find_name.search(f)
        return regmatch.groups()[0]

    rvt['sample'] = rvt.index.to_series().apply(samplify)
    
    rvt = rvt.groupby('sample').sum()
    rv = rvt.T.astype(int)
    lg.info("%d junctions observed", len(rv))
    
    # write to disk - raw count table
    lg.info('save to %s', args.datafile)
    rv.to_sql('counts', engine, if_exists='replace')


import logging
lg = logging.getLogger(__name__)

import warnings

import pandas as pd
from sqlalchemy import MetaData, create_engine
import sqlalchemy as sa
import leip
from path import Path

from nfj import util
    
@leip.flag('-c', '--counts')
@leip.arg('db')
@leip.command
def dbs(app, args):
    """
    List of databases
    """
    for tbl in Path(args.db).glob('*.pickle'):
        if args.counts:
            d = pd.read_pickle(tbl)
            print('%s\t%d' % (tbl.basename().replace('.pickle', ''), len(d)))
        else:
            print(tbl.basename().replace('.pickle', ''))

            
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('table')
@leip.command
def info(app, args):
    """
    Show table structure
    """
    from sqlalchemy import Index
    engine = create_engine(args.datafile)
    meta = MetaData()
    meta.reflect(engine)
    table = meta.tables[args.table]
    for col in table.columns:
        print('%s\t%s' % (col, col.type))
    for idx in table.indexes:
        print('IX %s\t%s' % (idx.name,
                          " ".join(map(str, idx.columns))))

@leip.arg('outfile', nargs='?')
@leip.arg('table')
@leip.arg('db')
@leip.command
def dump(app, args):
    d = util.load(args.db, args.table)
    if args.outfile:
        outfile = args.outfile
    else:
        outfile = args.table + ".tsv.gz"
    lg.info("Start writing %d records to %s", len(d), outfile) 
    d.to_csv(outfile, sep="\t")

    
@leip.arg('no', default=10, type=int, nargs='?')
@leip.arg('table')
@leip.arg('db')
@leip.command
def head(app, args):
    """
    Head on a table
    """
    d = util.load(args.db, args.table)
    print('# data shape:', d.shape)
    print(d.head(args.no))
    

@leip.arg('no', default=10, type=int, nargs='?')
@leip.arg('table')
@leip.arg('db')
@leip.command
def tail(app, args):
    """
    Head on a table
    """
    d = util.load(args.db, args.table)
    print('#data shape:', d.shape)
    print(d.tail(args.no))

        
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.command
def samples(app, args):
    """
    Show a list of samples
    """
    engine = create_engine(args.datafile)
    meta = MetaData()
    meta.reflect(engine)
    print("\t".join([str(x).replace('counts.', '')
                for x in meta.tables['counts'].columns
                if not x == 'counts.index']))


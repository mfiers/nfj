import warnings
warnings.filterwarnings("ignore")

import pandas as pd
from sqlalchemy import MetaData, create_engine
import sqlalchemy as sa
import leip

@leip.flag('-c', '--counts')
@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.command
def dbs(app, args):
    """
    List of databases
    """

    engine = create_engine(args.datafile)
    meta = MetaData()
    meta.reflect(engine)
    for tbl in meta.tables.keys():
        if args.counts:
            table = meta.tables[tbl]
            cnt = engine.execute(table.count()).first()[0]
            print('%s\t%s' % (tbl, cnt))
        else:
            print(tbl)

@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.command
def vacuum(app, args):
    """
    vacuum the database
    """
    engine = create_engine(args.datafile)
    engine.execute('VACUUM')

@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('table')
@leip.command
def drop(app, args):
    """
    drop a database
    """
    engine = create_engine(args.datafile)
    meta = MetaData()
    meta.reflect(engine)
    table = meta.tables[args.table]
    table.drop(engine)
    

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
        

@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('table')
@leip.command
def dump(app, args):
    engine = create_engine(args.datafile)
    d = pd.read_sql_table(args.table, engine)
    d.to_csv("%s.dump.tsv" % args.table, sep="\t")
                

@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('no', default=5, type=int, nargs='?')
@leip.arg('table')
@leip.command
def head(app, args):
    """
    Head on a table
    """
    engine = create_engine(args.datafile)
    meta = MetaData()
    meta.reflect(engine)
    table = meta.tables[args.table]
    head = engine.execute(table.select().limit(args.no))
    cols = table.columns
    cols = [str(x).replace('%s.' % args.table, '') for x in cols]
    print("\t".join(cols))
    for row in head:
        print("\t".join(map(str, row)))

        
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


@leip.arg('-d', '--datafile', default='sqlite:///nfj.db')
@leip.arg('no', default=5, type=int, nargs='?')
@leip.arg('table')
@leip.command
def tail(app, args):
    """
    Tail of a table
    """
    engine = create_engine(args.datafile)
    meta = MetaData()
    meta.reflect(engine)
    table = meta.tables[args.table]
    query = table.select().order_by(table.c.index.desc()).limit(args.no)
    head = engine.execute(query)
    cols = table.columns
    cols = [str(x).replace('%s.' % args.table, '') for x in cols]
    print("\t".join(cols))
    for row in head:
        print("\t".join(map(str, row)))
        


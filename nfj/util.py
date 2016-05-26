
import logging

import pandas as pd
from path import Path

lg = logging.getLogger(__name__)

def save(dbdir, **tables):
    dbdir = Path(dbdir).expanduser()
    dbdir.makedirs_p()
    lg.info("saving to db dir: %s" % dbdir)
    assert dbdir.isdir()
    for name, table in tables.items():
        lg.info("saving %s (%d records)", name, len(table))
        table.to_pickle(dbdir / '%s.pickle' % (name))
        table.describe().to_csv(
            dbdir / '%s.describe.tsv' % (name), sep="\t")

def load(db, name):
    dbdir = Path(db).expanduser()
    return pd.read_pickle(dbdir / ('%s.pickle' % name))
    
    

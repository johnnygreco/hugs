"""
Query hugs database and build sample. 
"""
from __future__ import division, print_function

from argparse import ArgumentParser
import pandas as pd
import hugs
from hugs.database.tables import Source
DB_FN = '/tigress/jgreco/hugs-dbs/20170709-173431/hsc-wide-patches-safe.db'

parser = ArgumentParser()
parser.add_argument('out_fn', type=str)
parser.add_argument('--db-fn', dest='db_fn', type=str, default=DB_FN)
parser.add_argument('--size-cut', dest='size_cut', type=float, default=2.5)
args = parser.parse_args()

engine = hugs.database.connect(args.db_fn)
session = hugs.database.Session()

num_sources = session.query(Source).count()
query = session.query(Source).filter(Source.flux_radius_i > args.size_cut)
num_big_sources = query.count()

frac_cut = (1.0 - num_big_sources/num_sources)*100
print('{} total source, {} big sources, {:.2f}% cut'.format(num_sources, 
    num_big_sources, frac_cut))

df = pd.read_sql(query.statement, engine)
df.to_csv(args.out_fn, index=False)

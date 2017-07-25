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
parser.add_argument('--size-cut-low', dest='size_cut_low', 
                    type=float, default=2.5)
parser.add_argument('--size-cut-high', dest='size_cut_high', 
                    type=float, default=20.0)
parser.add_argument('--do-color-cut', dest='do_color_cut', action='store_true')
args = parser.parse_args()

print('connecting to hugs database')
engine = hugs.database.connect(args.db_fn)
session = hugs.database.Session()

print('counting total sources')
num_sources = session.query(Source).count()

if args.do_color_cut:
    m, b = 0.7, 0.4
    color_line_lo =  lambda _x: m*_x - b
    color_line_hi =  lambda _x: m*_x + b
    gi = Source.mag_ap6_g - Source.mag_ap6_i - Source.A_g + Source.A_i
    gr = Source.mag_ap6_g - Source.mag_ap6_r - Source.A_g + Source.A_r
    query = session.query(Source)\
        .filter(Source.flux_radius_i > args.size_cut_low)\
        .filter(Source.flux_radius_i < args.size_cut_high)\
        .filter(gi > -0.1)\
        .filter(gi < 1.4)\
        .filter(color_line_lo(gi) < gr)\
        .filter(color_line_hi(gi) > gr)
else:
    query = session.query(Source)\
        .filter(Source.flux_radius_i > args.size_cut_low)\
        .filter(Source.flux_radius_i < args.size_cut_high)

print('converting query to pandas dataframe')
df = pd.read_sql(query.statement, engine)
num_big_sources = len(df)

frac_cut = (1.0 - num_big_sources/num_sources)*100
print('{} total sources, keeping {} sources, {:.2f}% cut'.format(num_sources, 
    num_big_sources, frac_cut))

print('writing catalog to file')
df.to_csv(args.out_fn, index=False)

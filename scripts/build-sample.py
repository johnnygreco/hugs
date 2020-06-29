"""
Query hugs database and build sample. 
"""
from __future__ import division, print_function

from argparse import ArgumentParser
import pandas as pd
import hugs
from hugs.database.tables import Source

parser = ArgumentParser()
parser.add_argument('out_fn', type=str)
parser.add_argument('--db-fn', dest='db_fn', type=str, required=True)
parser.add_argument('--size-cut-low', dest='size_cut_low', 
                    type=float, default=2.5)
parser.add_argument('--size-cut-high', dest='size_cut_high', 
                    type=float, default=20.0)
parser.add_argument('--do-color-cut', dest='do_color_cut', action='store_true')
parser.add_argument('--no-extinction', dest='no_ext', action='store_true')
args = parser.parse_args()

print('connecting to hugs database')
engine = hugs.database.connect(args.db_fn)
session = hugs.database.Session()

print('counting total sources')
num_sources = session.query(Source).count()

if args.do_color_cut:
    print('applying color and size cuts')
    m, b = 0.7, 0.4
    color_line_lo =  lambda _x: m*_x - b
    color_line_hi =  lambda _x: m*_x + b
    A_gi = 0.0 if args.no_ext else (Source.A_g - Source.A_i)
    A_gr = 0.0 if args.no_ext else (Source.A_g - Source.A_r)
    gi = Source.mag_ap6_g - Source.mag_ap6_i - A_gi
    gr = Source.mag_ap6_g - Source.mag_ap6_r - A_gr
    query = session.query(Source)\
        .filter(Source.flux_radius_50_i > args.size_cut_low)\
        .filter(Source.flux_radius_50_i < args.size_cut_high)\
        .filter(gi > -0.1)\
        .filter(gi < 1.4)\
        .filter(color_line_lo(gi) < gr)\
        .filter(color_line_hi(gi) > gr)
else:
    print('applying size cut')
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

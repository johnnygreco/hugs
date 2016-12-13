import pandas as pd
import hugs_pipe as hp
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('cat_fn', type=str)
parser.add_argument('-b', '--band', type=str, default='I')
parser.add_argument('-f', '--frame', type=int, help='ds9 frame', default=1)
args = hp.parse_args('patch', parser=parser)

cat = pd.read_csv(args.cat_fn)

if args.tract is None:
    assert args.patch is None
    tract, patch = cat.loc[0, ['tract', 'patch']]
else:
    tract, patch = args.tract, args.patch
data_id = {'tract':tract, 'patch':patch, 'filter': 'HSC-'+args.band.upper()}

viewer = hp.Viewer(data_id=data_id)
viewer.ds9_display_patch()
viewer.ds9_draw_ell(cat, frame=args.frame)

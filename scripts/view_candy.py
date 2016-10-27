from __future__ import division, print_function

from astropy.table import Table
import hugs_pipe

def main(cat_fn):
    cat = Table.read(cat_fn, format='ascii')
    hugs_pipe.viz.view_candy(cat)

if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser('view candidates')
    parser.add_argument('cat_fn', type=str, help='cat file name')
    args = parser.parse_args()
    main(args.cat_fn)

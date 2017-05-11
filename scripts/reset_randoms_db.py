from __future__ import print_function

import os
from argparse import ArgumentParser
from skyrandoms import SkyRandomsDatabase

parser = ArgumentParser()
parser.add_argument('db_fn', type=str)
args = parser.parse_args()

assert os.path.isfile(args.db_fn), 'db does not exist!'
db = SkyRandomsDatabase(args.db_fn)

print('setting all randoms to undetected in', args.db_fn)
db.set_all_undetected()

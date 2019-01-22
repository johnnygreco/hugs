import hugs
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-o', '--out-fn', dest='out_fn', required=True)
parser.add_argument('-c', '--config-fn', dest='config_fn', default=None, 
                    help='non-defulat synth parameters')
args = parser.parse_args()

if args.config_fn is not None:
    kwargs = hugs.utils.read_config(args.config_fn)
else: 
    # default parameters
    kwargs = dict(density=200,
                  ra_lim_list=[[0, 360]], 
                  dec_lim_list=[-15, 15], 
                  mu_range=[23, 28], 
                  r_eff_range=[3, 10], 
                  n_range=[0.3, 1.5], 
                  theta_range=[0, 180], 
                  min_sep=15, 
                  mu_type='average')

cat = hugs.synths.catalog.build_catalog_survey(**kwargs)
syncat = hugs.synths.catalog.SynthCat(catalog=cat)
syncat.write(args.out_fn)

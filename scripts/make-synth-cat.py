import hugs
from argparse import ArgumentParser

color_dict = dict(
    blues={'g_i': 0.47, 'g_r': 0.32},
    reds={'g_i': 0.82, 'g_r': 0.56},
    med={'g_i': 0.64, 'g_r': 0.43},
    random={'g_i': 'random', 'g_r': 'random'},
)

parser = ArgumentParser()
parser.add_argument('-o', '--out-fn', dest='out_fn', required=True)
parser.add_argument('--synth-color', required=True)
parser.add_argument('--q0', default=None, type=float)
parser.add_argument('-c', '--config-fn', dest='config_fn', default=None, 
                    help='non-defulat synth parameters')
args = parser.parse_args()

if args.config_fn is not None:
    kwargs = hugs.utils.read_config(args.config_fn)
else: 
    # default parameters
    kwargs = dict(density=200,
                  ra_lim_list=[[180, 200]], 
                  dec_lim_list=[[-5, 5]], 
                  mu_range=[23, 28], 
                  r_eff_range=[3, 10], 
                  n_range=[0.3, 1.5], 
                  theta_range=[0, 180], 
                  min_sep=15, 
                  mu_type='average')

kwargs.update(color_dict[args.synth_color])

if args.q0 is not None:
    kwargs['q0'] = args.q0

cat = hugs.synths.catalog.build_catalog_survey(**kwargs)
syncat = hugs.synths.catalog.GlobalSynthCat(catalog=cat)
syncat.write(args.out_fn)

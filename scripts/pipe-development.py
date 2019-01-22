import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
plt.style.use('../notebooks/jpg.mplstyle')
import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)

import hugs
from hugs.synths.catalog import SynthCat
from hugs.utils import project_dir
from hugs.pipeline import find_lsbgs
import tigerview


def _compare_param(param_in, param_out, param_name, type='against'):
    fig, ax = plt.subplots()

    if type == 'against':
        ax.plot(param_in, param_out, 'o')
        ax.plot([param_in.min(), param_in.max()], 
                [param_in.min(), param_out.max()], 'k-', lw=2)
        ax.set_xlabel(param_name + ' (injected)', fontsize=22)
        ax.set_ylabel(param_name + ' (recovered)', fontsize=22)
    elif type == 'diff':
        pass
    elif type == 'hist':
        ax.hist(param_in - param_out, bins='auto')
        ax.axvline(x=0, ls='--', lw=2, color='k')
        ax.set_xlabel(param_name, fontsize=22)

    ax.minorticks_on()
    plt.tight_layout()


def _plot_synths(recovered, missed):
    fig, ax = plt.subplots()

    ax.plot(recovered['r_e'], recovered['mu_0_g'], 'o')
    if len(missed) > 0:
        ax.plot(missed['r_e'], missed['mu_0_g'], 'kx')

    ax.set_xlabel(r'$r_\mathrm{eff}$')
    ax.set_ylabel(r'$\mu_0(g)$')

    ax.minorticks_on()
    plt.tight_layout()


def run(args):

    config = hugs.PipeConfig(tract=args.tract, patch=args.patch, 
                             config_fn=args.config_fn)

    results = find_lsbgs.run(config, reset_mask_planes=False)

    synth_cat_fn = hugs.utils.read_config(args.config_fn)['synth_cat_fn']
    synth_cat = SynthCat(cat_fn=synth_cat_fn).get_exp_synths(results.exp.i)

    synth_cat.rename_column('X0', 'x_image')
    synth_cat.rename_column('Y0', 'y_image')
    
    (match, match_synth), _  = hugs.cattools.xmatch(results.sources, 
                                                    synth_cat, 
                                                    max_sep=args.max_match_sep)
    source_match = results.sources[match]
    synth_match = synth_cat[match_synth]

    print('{} injected, {} recovered'.format(len(synth_cat), len(synth_match)))

    if not args.no_ds9:

        view = tigerview.ds9.Viewer()
        view.display_patch(exp=results.exp.i, frame=1, mask_trans=70)
        view.display_patch(exp=results.exp_clean, frame=2, mask_trans=70)

        if not args.no_ds9_sources:
            for src in source_match:
                view.display_source(src['ra'], src['dec'], 
                                    src['flux_radius_i']/0.168, 
                                    src['theta_image'], 
                                    src['ellipticity'])
            for src in synth_cat:
                view.display_source(src['ra'], src['dec'], src['r_e']/0.168, 
                                    src['theta'], src['ell'], color='red')


    _compare_param(source_match['mag_auto_i'], synth_match['m_i'], r'$m_i$')
    _compare_param(source_match['flux_radius_i'], synth_match['r_e'], 
                   r'$r_\mathrm{eff}$')
    _compare_param(source_match['mag_ap4_g'] - source_match['mag_ap4_i'], 
                   synth_match['g-i'], r'$\Delta(g-i)$', type='hist')
    _plot_synths(synth_match, synth_cat[~match_synth])


    # HACK: matplotlib crashing when scrolling or moving plot
    # https://github.com/matplotlib/matplotlib/issues/9637
    while True:
      try:
        plt.show()
      except UnicodeDecodeError:
        continue
      break


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('-t', '--tract', type=int, required=True)
    parser.add_argument('-p', '--patch', required=True)
    parser.add_argument('--no-ds9', dest='no_ds9', action='store_true')
    parser.add_argument('--no-ds9-sources', dest='no_ds9_sources', 
                        action='store_true')
    parser.add_argument('--min-match-sep', dest='max_match_sep', type=int,
                        help='maximum match separation in pixels', default=4)
    parser.add_argument(
        '-c', '--config-fn', dest='config_fn', 
        default=os.path.join(project_dir, 'pipe-configs/hugs-run-dev.yml'))
    args = parser.parse_args()

    run(args)

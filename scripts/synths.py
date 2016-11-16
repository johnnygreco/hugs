"""
Run hugs-pipe on data with synthetic UGDs.
"""
from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
from lsst.pipe.base import Struct
import hugs_pipe as hp

def calc_mask_bit_fracs(exp):
    mask = exp.getMaskedImage().getMask()
    msk_arr = mask.getArray()
    npix = float(msk_arr.size)
    getBitVal = mask.getPlaneBitMask

    npix_clean = (msk_arr & getBitVal('CLEANED') != 0).sum()
    npix_blend = (msk_arr & getBitVal('BLEND') != 0).sum()
    npix_bright = (msk_arr & getBitVal('BRIGHT_OBJECT') != 0).sum()

    return Struct(clean_frac=npix_clean/npix,
                  bright_frac=npix_bright/npix,
                  blend_frac=npix_blend/npix)

def main(config, pset_lims, outdir, num_synths=10, seed=None):

    synths_kwargs = {'num_synths': num_synths,
                     'seed': seed,
                     'pset_lims': pset_lims}
    tract, patch = config.data_id['tract'], config.data_id['patch']
    prefix = os.path.join(outdir, 'hugs-pipe-{}-{}'.format(tract, patch))

    results = hp.run(config,
                     debug_return=True,
                     inject_synths=True,
                     synths_kwargs=synths_kwargs)

    sources = results.sources
    sources.write(prefix+'.csv')

    sf = results.sf
    sf.write_cat(prefix+'-synths.csv')

if __name__=='__main__':
    args = hp.parse_args(os.path.join(hp.io, 'synth-results'))

    pset_lims = {'n': [0.6, 1.2],
                 'mu0_i': [23, 26],
                 'ell': [0.2, 0.2],
                 'r_e': [2, 5]}

    if args.group_id is None:
        tract, patch = args.tract, args.patch
        outdir = os.path.join(args.outdir, 'synths-{}-{}'.format(tract, patch))
        hp.utils.mkdir_if_needed(outdir)

        log_fn = os.path.join(outdir, 'hugs-pipe-synths.log')
        config = hp.Config(config_fn=args.config_fn, log_fn=log_fn)

        data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}
        config.set_data_id(data_id)
        main(config, pset_lims, outdir, 
             num_synths=args.num_synths, seed=args.seed)
    else:
        from astropy.table import Table
        regions_fn = 'cat_z0.065_Mh12.75-14.0_tracts_n_patches.npy'
        regions_fn = os.path.join(hp.io, regions_fn)
        regions_dict = np.load(regions_fn).item()
        regions = Table(regions_dict[args.group_id])

        outdir = args.group_dir
        hp.utils.mkdir_if_needed(outdir)
        log_fn = os.path.join(outdir, 'hugs-pipe-synths.log')
        config = hp.Config(log_fn=log_fn)

        print('searching in', len(regions), 'regions')
        for tract, patch in regions['tract', 'patch']:
            data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}
            config.set_data_id(data_id)
            main(config, pset_lims, outdir, 
                 num_synths=args.num_synths, seed=args.seed)
            config.reset_mask_planes()

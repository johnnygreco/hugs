"""
Run hugs-pipe on an HSC patch.
"""

import os
import numpy as np
import hugs_pipe
hugs_pipe_io = os.environ.get('HUGS_PIPE_IO')

def main(tract, patch, config, outdir):
    data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}
    hugs_pipe.utils.mkdir_if_needed(outdir)
    prefix = os.path.join(outdir, 'hugs-pipe-{}-{}'.format(tract, patch))
    if type(config)==str:
        config = hugs_pipe.Config(data_id=data_id, 
                                  config_fn=config, 
                                  log_fn=prefix+'.log')
    else:
        config.set_data_id(data_id)
    sources = hugs_pipe.run(config)
    sources.write(prefix+'.cat', format='ascii')

if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser('run hugs-pipe')
    parser.add_argument('-t', '--tract', type=int, help='HSC tract')
    parser.add_argument('-p', '--patch', type=str, help='HSC patch')
    parser.add_argument('-g', '--group_id', type=int, help='group id', 
                        default=None)
    parser.add_argument('-c', '--config', type=str, help='config file name', 
                        default=None)
    parser.add_argument('-o', '--outdir', type=str, help='output directory', 
                        default=hugs_pipe_io)
    args = parser.parse_args()
    if args.group_id is None:
        assert args.tract is not None
        assert args.patch is not None
        main(args.tract, args.patch, args.config, args.outdir)
    else:
        from astropy.table import Table
        outdir = os.path.join(args.outdir, 'group_'+str(args.group_id))
        log_fn = os.path.join(outdir, 'hugs-pipe.log')
        config = hugs_pipe.Config(log_fn=log_fn)
        regions_fn = 'cat_z0.065_Mh12.75-14.0_tracts_n_patches.npy'
        regions_fn = os.path.join(hugs_pipe_io, regions_fn)
        regions_dict = np.load(regions_fn).item()
        regions = Table(regions_dict[args.group_id])
        for tract, patch in regions['tract', 'patch']:
            main(tract, patch, config, outdir)

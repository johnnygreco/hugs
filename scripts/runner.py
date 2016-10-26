"""
Run hugs-pipe on an HSC patch.
"""

import os
import hugs_pipe

def main(tract, patch, config_fn, outdir):
    data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    prefix = os.path.join(outdir, 'hugs-pipe_{}_{}'.format(tract, patch))
    config = hugs_pipe.Config(data_id=data_id, 
                              config_fn=config_fn, 
                              log_fn=prefix+'.log')

    sources = hugs_pipe.run(config)
    sources.write(prefix+'.cat', format='ascii')

if __name__=='__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser('run hugs-pipe')
    parser.add_argument('tract', type=int, help='HSC tract')
    parser.add_argument('patch', type=str, help='HSC patch')
    parser.add_argument('-c', '--config', type=str, help='config file name', 
                        default=None)
    parser.add_argument('-o', '--outdir', type=str, help='output directory', 
                        default='/home/jgreco/hugs-pipe-out')
    args = parser.parse_args()
    main(args.tract, args.patch, args.config, args.outdir)

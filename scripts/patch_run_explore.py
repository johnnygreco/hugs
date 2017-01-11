import os
import hugs_pipe as hp

args = hp.parse_args(['patch', 'outdir', 'config_fn'])

assert (args.tract is not None) and (args.patch is not None)
tract, patch = args.tract, args.patch

outdir = os.path.join(args.outdir, 'testrun-{}-{}'.format(tract, patch))
hp.utils.mkdir_if_needed(outdir)

data_id = dict(tract=tract, patch=patch, filter='HSC-I')

prefix = os.path.join(outdir, 'hugs-pipe')

config = hp.Config(data_id=data_id, 
                   config_fn=args.config_fn, 
                   log_fn=prefix+'.log')

results = hp.run_use_sex(config, debug_return=True)

df = results.sources.to_pandas()

v1 = hp.Viewer(data_id=results.exposure)
disp1 = v1.create_ds9_display(frame=1)
disp1.setMaskPlaneColor('THRESH_HIGH', 'magenta')
disp1.setMaskPlaneColor('THRESH_LOW', 'yellow')
disp1.setMaskPlaneColor('CLEANED', 'white')
disp1.setMaskPlaneColor('SEX_SEG', 'blue')
v1.ds9_display_patch(frame=1, mask_trans=70, disp=disp1)
v1.ds9_draw_ell(df, frame=1, ellpars='sex')

v2 = hp.Viewer(data_id=results.exp_clean)
disp2 = v2.create_ds9_display(frame=2)
disp2.setMaskPlaneColor('THRESH_HIGH', 'magenta')
disp2.setMaskPlaneColor('THRESH_LOW', 'yellow')
disp2.setMaskPlaneColor('CLEANED', 'white')
disp2.setMaskPlaneColor('BLEND', 'black')
disp2.setMaskPlaneColor('SEX_SEG', 'blue')
v2.ds9_display_patch(frame=2, mask_trans=70, disp=disp2)
v2.ds9_draw_ell(df, frame=2, ellpars='sex')

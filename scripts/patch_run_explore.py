import os
import hugs_pipe as hp

args = hp.parse_args(['patch', 'outdir', 'config_fn'])

assert (args.tract is not None) and (args.patch is not None)
tract, patch = args.tract, args.patch

outdir = os.path.join(args.outdir, 'test-run-{}-{}'.format(tract, patch))

hp.utils.mkdir_if_needed(outdir)

prefix = os.path.join(outdir, 'hugs-pipe')

config = hp.Config(tract=tract,
                   patch=patch,
                   config_fn=args.config_fn, 
                   log_fn=prefix+'.log')

results = hp.run(config, debug_return=True)

fn = lambda name: os.path.join(outdir, name)

all_det = results.all_detections.to_pandas()
all_det.to_csv(fn('cat-all.csv'), index=False)

sources = results.sources.to_pandas()
sources.to_csv(fn('cat-verified.csv'))

candy = results.candy.to_pandas()
candy.to_csv(fn('cat-candy.csv'), index=False)

mask = results.exp_clean.getMaskedImage().getMask()
#mask.removeAndClearMaskPlane('SMOOTHED')
mask.removeAndClearMaskPlane('DETECTED')

v1 = hp.Viewer(data_id=results.exp['i'])
disp1 = v1.create_ds9_display(frame=1)
disp1.setMaskPlaneColor('THRESH_HIGH', 'magenta')
disp1.setMaskPlaneColor('THRESH_LOW', 'yellow')
disp1.setMaskPlaneColor('CLEANED', 'white')
disp1.setMaskPlaneColor('SEX_SEG', 'blue')
v1.ds9_display_patch(frame=1, mask_trans=70, disp=disp1)
v1.ds9_draw_ell(all_det, frame=1, ellpars='sex')
v1.ds9_draw_ell(sources, frame=1, ellpars='sex', scale=3, ec='red')
v1.ds9_draw_ell(candy, frame=1, ellpars='sex', scale=4, ec='green')

v2 = hp.Viewer(data_id=results.exp_clean)
disp2 = v2.create_ds9_display(frame=2)
disp2.setMaskPlaneColor('THRESH_HIGH', 'magenta')
disp2.setMaskPlaneColor('THRESH_LOW', 'yellow')
disp2.setMaskPlaneColor('CLEANED', 'white')
disp2.setMaskPlaneColor('BLEND', 'black')
disp2.setMaskPlaneColor('SEX_SEG', 'blue')
v2.ds9_display_patch(frame=2, mask_trans=70, disp=disp2)
v2.ds9_draw_ell(all_det, frame=2, ellpars='sex')
v2.ds9_draw_ell(sources, frame=2, ellpars='sex', scale=3, ec='red')
v2.ds9_draw_ell(candy, frame=2, ellpars='sex', scale=4, ec='green')

v3 = hp.Viewer(data_id=results.exp['g'])
disp3 = v2.create_ds9_display(frame=3)
disp3.setMaskPlaneColor('THRESH_HIGH', 'magenta')
disp3.setMaskPlaneColor('THRESH_LOW', 'yellow')
disp3.setMaskPlaneColor('CLEANED', 'white')
disp3.setMaskPlaneColor('BLEND', 'black')
disp3.setMaskPlaneColor('SEX_SEG', 'blue')
v3.ds9_display_patch(frame=3, mask_trans=70, disp=disp3)
v3.ds9_draw_ell(all_det, frame=3, ellpars='sex')
v3.ds9_draw_ell(sources, frame=3, ellpars='sex', scale=3, ec='red')
v3.ds9_draw_ell(candy, frame=3, ellpars='sex', scale=4, ec='green')

v4 = hp.Viewer(data_id=results.exp['r'])
disp4 = v2.create_ds9_display(frame=4)
disp4.setMaskPlaneColor('THRESH_HIGH', 'magenta')
disp4.setMaskPlaneColor('THRESH_LOW', 'yellow')
disp4.setMaskPlaneColor('CLEANED', 'white')
disp4.setMaskPlaneColor('BLEND', 'black')
disp4.setMaskPlaneColor('SEX_SEG', 'blue')
v4.ds9_display_patch(frame=4, mask_trans=70, disp=disp4)
v4.ds9_draw_ell(all_det, frame=4, ellpars='sex')
v4.ds9_draw_ell(sources, frame=4, ellpars='sex', scale=3, ec='red')
v4.ds9_draw_ell(candy, frame=4, ellpars='sex', scale=4, ec='green')

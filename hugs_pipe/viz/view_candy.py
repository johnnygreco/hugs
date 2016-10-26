"""
Display candidates with ds9. 
"""

from __future__ import division, print_function

import os
import lsst.afw.display as afwDisp
import lsst.afw.geom as afwGeom
from .. import imtools
hscdir = os.environ.get('HSC_DIR')

def view_candy(cat, butler=None):
    """
    View postage stamps of objects in catalog, which 
    must be the output of hugs-pipe.run.

    Parameters
    ----------
    cat : astropy.table.Table
        Object catalog.
    butler : lsst.daf.persistence.Butler, optional
        The butler. 
    """

    if butler is None:
        import lsst.daf.persistence
        butler = lsst.daf.persistence.Butler(hscdir)

    disp_1 = afwDisp.Display(1)
    disp_1.setMaskTransparency(100)
    disp_2 = afwDisp.Display(2)
    disp_2.setMaskTransparency(100)

    tract_old = None
    patch_old = None

    cat['size'] = 3*0.168*cat['semimajor_axis_sigma']
    cat['r_circ'] = 0.168*cat['equivalent_radius']

    for i, source in enumerate(cat): 

        print(cat['tract', 'patch', 'size', 'mu_3', 'mag_ell', 'r_circ'][i])

        tract = source['tract']
        patch = source['patch'][0]+','+source['patch'][-1]
        data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}

        if (tract_old!=tract) or (patch_old!=patch):
            exp = butler.get('deepCoadd_calexp', data_id, immediate=True)
            disp_1.erase()
            disp_1.mtv(exp)

        tract_old = tract
        patch_old = patch

        x, y = source['x_hsc'], source['y_hsc']
        cutout = imtools.get_cutout((x, y), 100, exp=exp)

        disp_2.erase()
        disp_2.mtv(cutout)
        a, b = source['semimajor_axis_sigma'], source['semiminor_axis_sigma']
        theta = source['orientation']

        with disp_2.Buffering():
            ell = afwGeom.ellipses.Axes(3*a, 3*b, theta)
            disp_2.dot(ell, x, y, ctype=afwDisp.CYAN)
            disp_2.dot('+', x, y, ctype=afwDisp.GREEN)   

        response = raw_input('\nPress enter to continue. Type stop to break: ')

        if response=='stop':
            break

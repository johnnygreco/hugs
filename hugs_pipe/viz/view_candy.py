"""
Display candidates with ds9. 
"""

from __future__ import division, print_function

import os
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
import lsst.afw.geom as afwGeom
from .. import imtools
from ..utils import pixscale
hscdir = os.environ.get('HSC_DIR')

__all__ = ['display_candies', 'view_stamps_on_tiger']


def display_candies(sources, data_id, frame=1):
    """
    """

    if type(data_id)==str:
        exp = afwImage.ExposureF(data_id)
    elif type(data_id)==afwImage.ExposureF:
        exp = data_id
    else:
        butler = lsst.daf.persistence.Butler(hscdir)
        exp = butler.get('deepCoadd_calexp', data_id, immediate=True)

    disp = afwDisp.Display(frame)
    disp.setMaskTransparency(75)
    disp.setMaskPlaneColor('THRESH_HIGH', 'magenta')
    disp.setMaskPlaneColor('THRESH_LOW', 'yellow')
    disp.setMaskPlaneColor('FAKE', 'green')
    disp.setMaskPlaneColor('CLEANED', 'white')
    disp.setMaskPlaneColor('BLEND', 'black')
    disp.mtv(exp)

    with disp.Buffering():
        for s in sources:
            x, y = s['x_hsc'], s['y_hsc']
            a, b = s['semimajor_axis_sigma'], s['semiminor_axis_sigma']
            theta = s['orientation']
            ell = afwGeom.ellipses.Axes(3*a, 3*b, theta)
            disp.dot(ell, x, y, ctype=afwDisp.CYAN)

    return disp


def view_stamps_on_tiger(cat, butler=None):
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

    for i, source in enumerate(cat): 

        print(cat['tract', 'patch', 'a_3_sig', 'mu_3', 'mag_ell', 'r_circ'][i])

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

"""
Display candidates with ds9. 
"""

from __future__ import division, print_function

import os
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
import lsst.afw.geom as afwGeom
import lsst.daf.persistence
from .. import imtools
from ..utils import pixscale
hscdir = os.environ.get('HSC_DIR')

__all__ = ['display_candies', 'display_cutout', 'view_stamps_on_tiger']


def display_candies(sources, data_id, frame=1, mask_trans=75):
    """
    Display image with candidates using ds9.
    """

    if type(data_id)==str:
        exp = afwImage.ExposureF(data_id)
    elif type(data_id)==afwImage.ExposureF:
        exp = data_id
    else:
        butler = lsst.daf.persistence.Butler(hscdir)
        exp = butler.get('deepCoadd_calexp', data_id, immediate=True)

    disp = afwDisp.Display(frame)
    disp.setMaskTransparency(mask_trans)
    disp.setMaskPlaneColor('THRESH_HIGH', 'magenta')
    disp.setMaskPlaneColor('THRESH_LOW', 'yellow')
    disp.setMaskPlaneColor('SYNTH', 'green')
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

def display_cutout(source, exp=None, data_id=None, frame=1, 
                   butler=None):
    """
    Display cutout of source with ds9.
    """

    if data_id:
        if butler is None:
            butler = lsst.daf.persistence.Butler(hscdir)
        exp = butler.get('deepCoadd_calexp', data_id, immediate=True)
    else:
        assert exp is not None

    disp = afwDisp.Display(frame)
    disp.setMaskTransparency(100)

    x, y = source['x_hsc'], source['y_hsc']    
    cutout = imtools.get_cutout((x, y), 100, exp=exp)
    disp.mtv(cutout)

    a, b = source['semimajor_axis_sigma'], source['semiminor_axis_sigma']
    theta = source['orientation']

    ell = afwGeom.ellipses.Axes(3*a, 3*b, theta)
    disp.dot(ell, x, y, ctype=afwDisp.CYAN)
    disp.dot('+', x, y, ctype=afwDisp.GREEN)

    return disp


def view_stamps_on_tiger(cat, butler=None, frame=1):
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
        butler = lsst.daf.persistence.Butler(hscdir)

    disp = afwDisp.Display(frame)
    disp.setMaskTransparency(100)

    tract_old = None
    patch_old = None

    for i, source in enumerate(cat): 

        print(cat['tract', 'patch', 'a_3_sig', 'mu_3_i', 'mag_ell_i',
                  'r_circ_ell'][i])

        tract = source['tract']
        patch = source['patch'][0]+','+source['patch'][-1]
        data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}

        if (tract_old!=tract) or (patch_old!=patch):
            exp = butler.get('deepCoadd_calexp', data_id, immediate=True)

        tract_old = tract
        patch_old = patch

        x, y = source['x_hsc'], source['y_hsc']
        cutout = imtools.get_cutout((x, y), 100, exp=exp)

        disp.erase()
        disp.mtv(cutout)
        a, b = source['semimajor_axis_sigma'], source['semiminor_axis_sigma']
        theta = source['orientation']

        with disp.Buffering():
            ell = afwGeom.ellipses.Axes(3*a, 3*b, theta)
            disp.dot(ell, x, y, ctype=afwDisp.CYAN)
            disp.dot('+', x, y, ctype=afwDisp.GREEN)   

        response = raw_input('\nPress enter to continue. Type stop to break: ')

        if response=='stop':
            break

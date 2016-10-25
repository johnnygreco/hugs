from __future__ import division, print_function

import os
import lsst.daf.persistence
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
import lsst.afw.geom as afwGeom
hscdir = os.environ.get('HSC_DIR')

def get_cutout(center, size, exp=None, data_id=None, butler=None):
    """
    Generate a cutout of exposure. Most give exposure object or 
    data_id and, optionally, a butler.
     
    Parameters
    ----------
    center : tuple
        Center of desired cutout in the tract
        and patch system.
    size : float
        Size to grow bbox in all directions.
    exp : lsst.afw.image.ExposureF, optional
        Exposure from which to get cutout.    
    data_id : dict, optional
        HSC data ID.
    butler : lsst.daf.persistence.Butler, optional
        the butler.

    Returns
    -------
    cutout : lsst.afw.image.ExposureF
        Desired exposure cutout.
    """

    # get exposure
    if exp is None:
        if butler is None:
            butler = lsst.daf.persistence.Butler(hscdir)
        exp = butler.get('deepCoadd_calexp', data_id, immediate=True)

    # generate bbox and get cutout
    center = afwGeom.Point2I(int(center[0]), int(center[1]))
    bbox = afwGeom.Box2I(center, center)
    bbox.grow(size)
    bbox.clip(exp.getBBox())
    cutout = exp.Factory(exp, bbox, afwImage.PARENT)

    return cutout

from __future__ import division, print_function

import numpy as np
import lsst.afw.geom
import lsstutils
from ..log import logger


DEFAULT_PLANES = ['BRIGHT_OBJECT', 'NO_DATA']


def get_mask_array(exp, planes=DEFAULT_PLANES):
    mask = exp.getMaskedImage().getMask()
    arr = np.zeros(mask.getArray().shape, dtype=bool)
    for p in planes:
        if p in mask.getMaskPlaneDict().keys():
            arr |= mask.getArray() & mask.getPlaneBitMask(p) != 0
    return arr


def find_masked_synths(synth_cat, exp):

    afwcoords = lsstutils.make_afw_coords(synth_cat['ra', 'dec'])

    # get mask array and find masked synths
    xy0 = exp.getXY0()
    wcs = exp.getWcs()
    bbox = exp.getBBox()
    mask_arr = get_mask_array(exp)
    masked= []
    for coord in afwcoords:
        pixel = wcs.skyToPixel(coord)
        if bbox.contains(lsst.afw.geom.Point2I(pixel)):
            j, i = pixel - xy0
            masked.append(int(mask_arr[int(i), int(j)]))
        else:
            logger.warn('synth not in exposure')
            masked.append(0)

    masked = np.array(masked)

    msg = '{} out of {} synths were masked'
    logger.info(msg.format(np.sum(masked), len(synth_cat)))

    return masked

from __future__ import division, print_function

import numpy as np
import lsstutils
from skyrandoms import SkyRandomsDatabase

__all__ = ['get_mask_array', 'find_randoms_in_footprint']
DEFAULT_PLANES = ['CLEANED', 'BRIGHT_OBJECT', 'BLEND', 'NO_DATA']

def get_mask_array(exp, planes=DEFAULT_PLANES):
    mask = exp.getMaskedImage().getMask()
    arr = np.zeros(mask.getArray().shape, dtype=bool)
    for p in planes: 
        if p in mask.getMaskPlaneDict().keys():
            arr |= mask.getArray() & mask.getPlaneBitMask(p) != 0
    return arr 


def find_randoms_in_footprint(db_fn, exp, return_db=True):
    """
    Find random points that fall within patch footprint and 
    are not masked by hugs-pipe. 

    Parameters
    ----------
    db_fn : string
        Database file name.
    exp : lsst.afw.ExposureF
        Exposure with WCS and mask.
    update_db : bool
        If True, update database. Else, return dataframe

    Returns
    -------
    df : pandas.DataFrame
        SkyRandoms table with updated detected column.
    db : SkyRandomsDatabase, it return_db=True
        Database manager. 
    """

    # load randoms database
    db = SkyRandomsDatabase(db_fn)

    # query database 
    ra_lim, dec_lim = lsstutils.get_exp_sky_limits(exp)
    df = db.query_region(ra_lim, dec_lim)
    afwcoords = lsstutils.make_afw_coords(df[['ra', 'dec']].values)

    # get mask array and find detected randoms
    xy0 = exp.getXY0()
    wcs = exp.getWcs()
    mask_arr = get_mask_array(exp)
    detected = []
    for coord in afwcoords:
        j, i = wcs.skyToPixel(coord) - xy0
        detected.append(int(not mask_arr[int(j), int(i)]))
    df['detected'] = detected

    return (df, db) if return_db else df

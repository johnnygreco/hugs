from __future__ import division, print_function

import numpy as np
from scipy.spatial import cKDTree
from astropy.table import Table
from .utils import check_random_state, check_kwargs_defaults

__all__ = ['find_duplicates', 
           'get_random_subset',
           'remove_duplicates', 
           'xmatch']

def find_duplicates(cat, max_sep=0.7, ra_col='ra', dec_col='dec'):
    """
    Find duplicate sources in catalog.

    Parameters
    ----------
    cat : pandas DataFrame
        Source catalog. 
    max_sep : float, optional
        Max separation to be considered the same source (in arcsec).
    ra_col, dec_col : string, optional
        The name of the ra and dec columns, which should be in degrees. 

    Returns
    -------
    ind : ndarray
        Index of matching pairs (i,j) with i<j. 
    """
    kdt = cKDTree(cat[[ra_col, dec_col]].values)
    ind = kdt.query_pairs(max_sep/3600.0, output_type='ndarray')
    return ind


def get_random_subset(cat, size, random_state=None):
    """
    Get random subset of source catalog.

    Parameters
    ----------
    cat : pandas DataFrame
        Source catalog.
    size : int
        Size of the random sample. 
    random_state : int, RandomState instance or None, optional 
        If int, random_state is the rng seed.
        If RandomState instance, random_state is the rng.
        If None, the rng is the RandomState instance used by np.random.

    Returns
    -------
    subset : pandas DataFrame
        The random subset.
    """
    rng = check_random_state(random_state)
    random_rows = rng.choice(len(cat), size, replace=False)
    subset = cat.iloc[random_rows]
    return subset


def remove_duplicates(cat, inplace=True, **kwargs):
    """
    Remove duplicate entries from in source catalog.

    Parameters
    ----------
    cat : pandas DataFrame
        Source catalog. 
    **kwargs : dict
        Keywords for find_duplicates.
    """
    ind = find_duplicates(cat, **kwargs)
    return cat.drop(cat.index[ind[:,1]], inplace=inplace)


def xmatch(cat_1, cat_2, xy_cols=['X_IMAGE', 'Y_IMAGE'], max_sep=5):
    """
    Crossmatch catalogs.
    """
    if type(cat_1)==Table:
        coords_1 = cat_1[xy_cols].to_pandas().values
    else:
        coords_1 = cat_1[xy_cols].values
    if type(cat_2)==Table:
        coords_2 = cat_2[xy_cols].to_pandas().values
    else:
        coords_2 = cat_2[xy_cols].values
    kdt = cKDTree(coords_1)
    dist, idx = kdt.query(coords_2)
    match_2 = dist < max_sep
    match_1 = idx[match_2]
    mismatch_2 = ~match_2
    mismatch_1 = idx[mismatch_2]
    match_masks = (match_1, match_2)
    mismatch_masks = (mismatch_1, mismatch_2)
    return match_masks, mismatch_masks 

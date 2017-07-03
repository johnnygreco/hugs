from __future__ import division, print_function

import numpy as np
from scipy.spatial import cKDTree
from astropy.table import Table
from .utils import check_random_state, check_kwargs_defaults

__all__ = ['cutter', 
           'find_duplicates', 
           'get_random_subset',
           'max_r_vir_mask',
           'remove_duplicates', 
           'xmatch']

##################################
# Default selection cuts
##################################

MIN_CUTS = {'FWHM_IMAGE': 25} 
MAX_CUTS = {'num_edge_pix': 1}

def cutter(cat, min_cuts=MIN_CUTS, max_cuts=MAX_CUTS, 
           inplace=False, group_id=None, max_r_vir=2.0, 
           cut_duplicates=False, logger=None, **kwargs):
    """
    Make selection cuts on input catalog.

    Parameters
    ----------
    cat : pandas DataFrame
        Input catalog.
    min_cuts : dict, optional
        Minimum params and values to cut.
    max_cuts : dict, optional
        Maximum params and values to cut.
    inplace: bool, optional
        If True, cut catalog inplace.
    group_id : int, optional
        Galaxy group id. If given, will impose max 
        virial radius cut.
    max_r_vir : float, optional
        Max number of virial radii to be considered a candidate. If None, 
        no cut will be made.
    cut_duplicates : bool
        If True, find and remove duplicate entries using remove_duplicates, 
        which can be tuned with kwargs.

    Returns 
    -------
    cut_cat : ndarray or None
        Cut catalog. Will be None if inplace=True
    """

    if logger:
        log = logger.info
    else:
        log = print

    log('cuts: {} objects in cat before cuts'.format(len(cat)))

    min_mask = np.ones(len(cat), dtype=bool)
    for key, min_val in min_cuts.items():
        if min_val is not None:
            cut = cat[key].values <= min_val
            min_mask[cut] = False
            log('cuts: will cut {} objects with {} <= {}'.format(
                (cut).sum(), key, min_val))

    max_mask = np.ones(len(cat), dtype=bool)
    for key, max_val in max_cuts.items():
        if max_val is not None:
            cut = cat[key].values >= max_val
            max_mask[cut] = False
            log('cuts: will cut {} objects with {} >= {}'.format(
                (cut).sum(), key, max_val))

    mask = min_mask & max_mask

    if group_id and max_r_vir:
        cut = max_r_vir_mask(cat, group_id, max_r_vir)
        mask &= cut
        log('cuts: will cut {} objects outside {} r_vir'.format(
            (cut).sum(), max_r_vir))

    log('cuts: {} objects in cat after cuts'.format(mask.sum()))

    cut_cat = cat.drop(cat.index[~mask], inplace=inplace)

    if cut_duplicates:
        if inplace:
            cut_cat = remove_duplicates(cat, inplace=inplace)
            log('cuts: {} objects after removing duplicates'.format(len(cat)))
        else:
            cut_cat = remove_duplicates(cut_cat, inplace=inplace)
            log('cuts: {} objects after removing duplicates'.format(
                len(cut_cat)))

    return cut_cat


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


def max_r_vir_mask(cat, group_id, max_r_vir=2.0):
    """
    If we have the group_id, we can cut objects 
    that are greater than max_r_vir from the group. 
    """
    from hugs.datasets import yang
    from toolbox.astro import angsep
    from toolbox.cosmo import Cosmology
    cosmo = Cosmology()
    props = ['ra', 'dec', 'z', 'Mh_Lest']
    ra, dec, z, logMh = yang.get_group_prop(group_id, props)
    D_A = cosmo.D_A(z)
    r180 = yang.r180(logMh, z)
    theta_180 = (r180/D_A)*180.0/np.pi
    seps = angsep(ra, dec, cat['ra'], cat['dec'], sepunits='degree')
    mask = seps < max_r_vir*theta_180
    return mask


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

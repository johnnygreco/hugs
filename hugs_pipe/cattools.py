from __future__ import division, print_function

import numpy as np
from .utils import check_random_state, check_kwargs_defaults

__all__ = ['cutter', 
           'find_duplicates', 
           'get_random_subset',
           'max_r_vir_mask',
           'remove_duplicates']

##################################
# Default selection cuts
##################################

MIN_CUTS = {'a_2_sig': 2.0,
            'mu_2_i': 23.0}
MAX_CUTS = {'num_edge_pix': 1}


def cutter(cat, min_cuts=MIN_CUTS, max_cuts=MAX_CUTS, verbose=True, 
           inplace=False, group_id=None, max_r_vir=2.0, 
           cut_duplicates=False, **kwargs):
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
    verbose : bool, optional
        If True, print lots of info.
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

    if verbose:
        print(len(cat), 'objects in cat before cuts')

    min_mask = np.ones(len(cat), dtype=bool)
    for key, min_val in min_cuts.items():
        if min_val is not None:
            cut = cat[key].values <= min_val
            min_mask[cut] = False
            if verbose:
                print('will cut {} objects with {} <= {}'.format(
                    (cut).sum(), key, min_val))

    max_mask = np.ones(len(cat), dtype=bool)
    for key, max_val in max_cuts.items():
        if max_val is not None:
            cut = cat[key].values >= max_val
            max_mask[cut] = False
            if verbose:
                print('will cut {} objects with {} >= {}'.format(
                    (cut).sum(), key, max_val))

    mask = min_mask & max_mask

    if group_id and max_r_vir:
        cut = max_r_vir_mask(cat, group_id, max_r_vir)
        mask &= cut
        if verbose:
            print('will cut {} objects outside {} r_vir'.format(
                  (cut).sum(), max_r_vir))

    if verbose:
        print(mask.sum(), 'objects in cat after cuts')

    cut_cat = cat.drop(cat.index[~mask], inplace=inplace)

    if cut_duplicates:
        if inplace:
            cut_cat = remove_duplicates(cat, inplace=inplace)
            if verbose:
                print(len(cat), 'objects after removing duplicates')
        else:
            cut_cat = remove_duplicates(cut_cat, inplace=inplace)
            if verbose:
                print(len(cut_cat), 'objects after removing duplicates')

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
    from scipy.spatial import cKDTree
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


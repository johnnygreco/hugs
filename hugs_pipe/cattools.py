from __future__ import division, print_function

import numpy as np

__all__ = ['find_duplicates', 'get_random_subset', 'remove_duplicates']


def find_duplicates(cat, max_sep=0.7, ra_col='ra', dec_col='dec'):
    """
    Find duplicate sources in catalog.

    Parameters
    ----------
    cat : pandas DataFrame
        Source catalog. 
    max_sep : float
        Max separation to be considered the same source (in arcsec).
    ra_col, dec_col : string
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


def get_random_subset(cat, size):
    """
    Get random subset of source catalog.

    Parameters
    ----------
    cat : pandas DataFrame
        Source catalog.
    size : int
        Size of the random sample. 

    Returns
    -------
    subset : pandas DataFrame
        The random subset.
    """
    random_rows = np.random.choice(len(cat), size, replace=False)
    subset = cat.iloc[random_rows]
    return subset


def remove_duplicates(cat, **kwargs):
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
    cat.drop(cat.index[ind[:,1]], inplace=True)

import os
import numpy as np
from scipy.special import gammaincinv, gamma
from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord, concatenate
import lsst.afw.geom
from ..log import logger
from ..utils import check_random_state, angsep
from ..utils import ra_dec_to_xyz, angular_dist_to_euclidean_dist 
from ..utils import euclidean_dist_to_angular_dist, read_config, solid_angle


__all__ = ['random_radec', 'synthetic_sersics', 'build_catalog_field',
           'build_catalog_survey', 'random_positions', 'generate_patch_cat', 
           'GlobalSynthCat']


def random_radec(nsynths, ra_lim=[0, 360], dec_lim=[-90, 90],
                 random_state=None, **kwargs):
    """
    Generate random ra and dec points within a specified range.
    All angles in degrees.
    Parameters
    ----------
    nsynths : int
        Number of random points to generate.
    ra_lim : list-like, optional
        ra limits.
    dec_lim : list-like, optional
        dec limits.
    random_state : `None`, int, list of ints, or `numpy.random.RandomState`
        If ``seed`` is `None`, return the `~numpy.random.RandomState`
        singleton used by ``numpy.random``.  If ``seed`` is an `int`,
        return a new `~numpy.random.RandomState` instance seeded with
        ``seed``.  If ``seed`` is already a `~numpy.random.RandomState`,
        return it.  Otherwise raise ``ValueError``.

    Returns
    -------
    points : 2d ndarray 
        Random ra and dec points in degrees.
    """
    rng = check_random_state(random_state)
    ra_lim = np.deg2rad(np.asarray(ra_lim))
    dec_lim = np.deg2rad(np.asarray(dec_lim))

    zlim = np.sin(dec_lim)
    z = zlim[0] + zlim.ptp() * rng.uniform(size=int(nsynths))
    ra = ra_lim[0] + ra_lim.ptp() * rng.uniform(size=int(nsynths))
    dec = np.arcsin(z)
    ra, dec = np.rad2deg(ra), np.rad2deg(dec)
    points = np.array([ra, dec]).T

    return points


def random_positions(image_shape=(501, 501), nsynths=10, edge_buffer=20,
                     random_state=None, **kwargs):

     rng = check_random_state(random_state)
     ydim, xdim = image_shape[0] - edge_buffer, image_shape[1] - edge_buffer
     x = rng.randint(edge_buffer, image_shape[1] - edge_buffer, nsynths)
     y = rng.randint(edge_buffer, image_shape[0] - edge_buffer, nsynths)

     return x, y


def synthetic_sersics(mu_range=[23, 28], r_eff_range=[3, 15], 
                      n_range=[0.3, 1.5], ell_range=[0., 0.65], 
                      theta_range=[0, 180], nsynths=100, random_state=None,  
                      master_band='g', mu_type="central", 
                      g_i=0.6, g_r=0.4, **kwargs):
    """
    Generate catalog of Sersic galaxies.
    

    Notes
    -----
    full sample:
    median g-i = 0.64
    median g-r = 0.43

    reds:
    median g-i = 0.82
    median g-r = 0.56

    blues:
    median g-i = 0.47
    median g-r = 0.32
    """

    size = int(nsynths)

    rng = check_random_state(random_state)
    r_eff = rng.uniform(*r_eff_range, size=size)    
    sersic_n = rng.uniform(*n_range, size=size)    
    ell = rng.uniform(*ell_range, size=size)    
    theta = rng.uniform(*theta_range, size=size)    
    
    b_a = 1 - ell
    b_n = gammaincinv(2.*sersic_n, 0.5)
    f_n = gamma(2*sersic_n)*sersic_n*np.exp(b_n)/b_n**(2*sersic_n)

    mu = rng.uniform(*mu_range, size=size)    

    if mu_type=='central':
        mu_0 = mu
        mu_e = mu_0 + 2.5*b_n/np.log(10)
        mu_e_ave = mu_e - 2.5*np.log10(f_n)
    elif mu_type=='average':
        mu_e_ave = mu
        mu_e = mu_e_ave + 2.5*np.log10(f_n)
        mu_0 = mu_e - 2.5*b_n/np.log(10) 
    else:
        raise Exception(mu_type+' is not a valid mu type')

    r_circ = r_eff*np.sqrt(b_a)
    A_eff = np.pi*r_circ**2
    m_tot = mu_e_ave - 2.5*np.log10(2*A_eff)    

    cat = {'m_' + master_band: m_tot, 
            'mu_0_' + master_band: mu_0, 
            'mu_e_ave_' + master_band: mu_e_ave}

    if master_band == 'g':
        # write i band
        cat['m_i'] = m_tot - g_i
        cat['mu_0_i'] = mu_0 - g_i
        cat['mu_e_ave_i'] = mu_e_ave - g_i
        # write r band
        cat['m_r'] = m_tot - g_r
        cat['mu_0_r'] = mu_0 - g_r
        cat['mu_e_ave_r'] = mu_e_ave - g_r
    elif master_band == 'i':
        # write g band
        cat['m_g'] = m_tot + g_i
        cat['mu_0_g'] = mu_0 + g_i
        cat['mu_e_ave_g'] = mu_e_ave + g_i
        # write r band
        cat['m_r'] = m_tot + g_i - g_r
        cat['mu_0_r'] = mu_0 + g_i - g_r
        cat['mu_e_ave_r'] = mu_e_ave + g_i - g_r
    else:
        raise Exception('master_band must be g or i')

    cat = Table(cat)
   
    cat['theta'] = theta
    cat['PA'] = theta - 90
    cat['r_e'] = r_eff
    cat['ell'] = 1 - b_a 
    cat['n'] = sersic_n
    cat['g-i'] = g_i
    cat['g-r'] = g_r

    return cat


def build_catalog_field(min_sep=None, **kwargs):
    """
    Build a catalog of synthetic sersic galaxies within a single field 
    (which may be as big as the entire sky).
    """
    nsynths = kwargs.pop('nsynths')

    if min_sep is None:
        coords = random_radec(nsynths=nsynths, **kwargs)
    else:
        count = 0
        while count < nsynths:
            ra, dec = random_radec(nsynths=1, **kwargs)[0]
            if count == 0:
                coords = np.array([[ra, dec]])
                count += 1
            else:
                sep = angsep(ra, dec, coords[:, 0], coords[:, 1])
                if (sep < min_sep).sum() == 0:
                    coords = np.append(coords, [[ra, dec]], axis=0)
                    count += 1

    catalog = synthetic_sersics(nsynths=nsynths, **kwargs)
    catalog['ra'] = coords[:, 0]
    catalog['dec'] = coords[:, 1]

    return catalog


def build_catalog_survey(density, ra_lim_list, dec_lim_list, mu_range, 
                         r_eff_range, n_range, theta_range, min_sep, **kwargs):
    """
    Build a catalog of synthetic sersic galaxies for a survey that is specified
    by a list of ra and dec limits. 

    Parameters
    ----------
    density : int
        Number density of sources per square degree
    ra_lim_list : list 
        List of ra limits (e.g., [[ra_1_min, ra_1_max], [ra_2_min, ra_2_max]]) 
    dec_lim_list : list 
        List of dec limits
    mu_range : list
        Range of surface brightness values for mocks
    r_eff_range : list
        Range of effective radii in arcsec
    n_range : list
        Range of sersic n values
    theta_range : list
        Range of position angles
    min_sep : float
        Minimum separation between mocks in arcsec
    kwargs : 
        Optional keywords for synthetic_sersics and random_radec

    Returns
    -------
    cat : astropy.Table
        The mock catalog
    """
    cat = []
    for ra_lim, dec_lim in zip(ra_lim_list, dec_lim_list):
        area = solid_angle(ra_lim, dec_lim)
        logger.info('generating catalog over {:.2f} deg^2'.format(area))
        nsynths = int(np.ceil(density * area))
        logger.info('nsynths = {}'.format(nsynths))
        kws = dict(
            nsynths=nsynths, ra_lim=ra_lim, dec_lim=dec_lim, mu_range=mu_range,
            r_eff_range=r_eff_range, n_range=n_range, theta_range=theta_range, 
            min_sep=min_sep
        )
        kws.update(kwargs)
        cat.append(build_catalog_field(**kws))
    cat = vstack(cat)
    return cat


class GlobalSynthCat(object):
    """
    A class for synthetic catalogs with a KDTree attribute to allow 
    for super fast queries. 
    """
    
    def __init__(self, cat_fn=None, catalog=None, cat_params={'nsynths':100}):
        from sklearn.neighbors import KDTree

        if cat_fn is not None:
            self.cat = Table.read(cat_fn)
        elif catalog is not None:
            self.cat = catalog
        else:
            self.cat = build_synthetic_galaxy_catalog(**cat_params)

        self.cat['synth_id'] = np.arange(1, len(self.cat) + 1)

        xyz = ra_dec_to_xyz(self.cat['ra'], self.cat['dec'])
        self.kdt = KDTree(np.asarray(xyz).T)

    def query_radius(self, ra, dec, r):
        """
        Search for sources around coordinate within circle of 
        radius r in arcseconds. 
        """
        xyz = np.array(ra_dec_to_xyz(ra, dec)).T.reshape(1, -1)
        idx = self.kdt.query_radius(
            xyz, angular_dist_to_euclidean_dist(r / 3600.0), 
            count_only=False, return_distance=False)[0]
        return self.cat[idx]

    def get_exp_synths(self, exp, search_radius=720):
        """
        Get synths in that fall within the given exposure.
        """

        wcs = exp.getWcs()
        xc, yc =  exp.getDimensions()//2 +  exp.getXY0()
        coord = wcs.pixelToSky(lsst.afw.geom.Point2D(xc, yc))
        ra_c, dec_c = coord.getRa().asDegrees(), coord.getDec().asDegrees()
        cat = self.query_radius(ra_c, dec_c, search_radius).copy()

        if len(cat) > 0:
            mask = np.zeros(len(cat), dtype=bool)
            cat['x'] = -1
            cat['y'] = -1
            
            for i, src in enumerate(cat):
                sky_coord = lsst.afw.geom.SpherePoint(
                    src['ra'] * lsst.afw.geom.degrees, 
                    src['dec'] * lsst.afw.geom.degrees)
                xy_coord = wcs.skyToPixel(sky_coord)
                if exp.getBBox().contains(lsst.afw.geom.Point2I(xy_coord)):
                    mask[i] = True
                    x0, y0 = xy_coord - exp.getXY0()
                    cat[i]['x'] = x0
                    cat[i]['y'] = y0

            cat = cat[mask]
            self.set_injected(cat['synth_id'])
        
        return cat

    def set_injected(self, synth_ids):
        idx = np.asarray(synth_ids) - 1
        self.cat['injected'][idx] = True

    def write(self, fn):
        self.cat.write(fn, overwrite=True)


def generate_patch_cat(nsynths, image_shape, edge_buffer, sersic_params={}, 
                       random_state=None, min_pixel_sep=150, wcs=None):

    x, y = random_positions(image_shape, 1, edge_buffer, random_state)

    cat_size = 1
    while cat_size < nsynths:
        _x, _y = random_positions(image_shape, 1, edge_buffer, random_state)
        seps = np.sqrt((_x - x)**2 + (_y - y)**2)
        if (seps < min_pixel_sep).sum() == 0:
            x = np.concatenate([x, _x])
            y = np.concatenate([y, _y])
            cat_size += 1

    sersic_params['nsynths'] = nsynths
    catalog = synthetic_sersics(**sersic_params)
    catalog['x'] = x
    catalog['y'] = y

    return catalog

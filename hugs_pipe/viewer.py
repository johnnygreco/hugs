"""
View HSC patches and hugs-pipe sources with ds9.
"""

from __future__ import division, print_function

import os
import matplotlib.pyplot as plt
import lsst.pipe.base
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
import lsst.afw.display.rgb as afwRgb
import lsst.afw.geom as afwGeom
from . import imtools
from .utils import pixscale, get_exposure, check_astropy_to_pandas
DATADIR = os.environ.get('HSC_DIR')

__all__ = ['Viewer']

class Viewer(object):

    def __init__(self, data_id, butler=None, 
                 data_dir=DATADIR, hsc_rgb_bands='IRG'):
        """
        Parameters
        ----------
        data_id : str, afwImage.ExposureF, or dict
            Exposure info. Can be file name, an exposure object, 
            or dict of HSC data ID.
        butler : lsst.daf.persistence.Butler, optional
            The Butler.
        data_dir : str, optional
            Location of exposure data.
        hsc_rgb_bands : str, optional
            HSC bands in 'RGB' order for visualization.
        """

        self._data_id = data_id
        self._butler = butler
        self.data_dir = data_dir
        self.frames = {}
        self.exp = get_exposure(data_id, self.butler, data_dir)
        self._rgb_images = None
        self.hsc_rgb_bands = hsc_rgb_bands

    @property
    def butler(self):
        """
        The Butler.
        """
        if self._butler is None:
            import lsst.daf.persistence
            self._butler = lsst.daf.persistence.Butler(self.data_dir)
        return self._butler

    @property
    def rgb_images(self):
        """
        Color data for visualization.

        Parameters
        ----------
        hsc_bands : str
            HSC filters in 'RGB' order.  
        reload : bool
            If true, reload if already loaded. 
        """

        if self._rgb_images is None:
            self._rgb_images = {}
            for band in self.hsc_rgb_bands:
                _id = self.get_data_id()
                _id['filter'] = 'HSC-'+band
                exp = self.butler.get('deepCoadd_calexp', _id, immediate=True)
                self._rgb_images.update({band: exp.getMaskedImage()})
        return self._rgb_images

    def get_data_id(self):
        return self._data_id.copy()

    def ds9_display_patch(self, frame=1, mask_trans=100):
        """
        Display patch with candidates using ds9.

        Parameters
        ----------
        frame : int, optional
            ds9 frame.
        mask_trans : float, optional
            Mask transparency
        """

        disp = afwDisp.Display(frame)
        disp.mtv(self.exp)
        disp.setMaskTransparency(mask_trans)
        frame_props = lsst.pipe.base.Struct(
            disp=disp, bbox=self.exp.getBBox(afwImage.PARENT)) 
        self.frames.update({frame:frame_props})

    def ds9_display_cutout(self, coord_hsc, frame=1, mask_trans=100, size=100):
        """
        Display cutout of source with ds9.

        Parameters ----------
        coord_hsc : tuple
            Central coordinate in HSC tract system.
        frame : int, optional
            ds9 frame.
        mask_trans : float, optional
            Mask transparency
	size : int, optional
             Size to grow bbox in all directions.
        """

        disp = afwDisp.Display(frame)
        cutout = imtools.get_cutout(coord_hsc, size, exp=self.exp)
        disp.mtv(cutout)
        disp.setMaskTransparency(mask_trans)
        frame_props = lsst.pipe.base.Struct(
            disp=disp, bbox=cutout.getBBox(afwImage.PARENT)) 
        self.frames.update({frame:frame_props})

    def ds9_draw_ell(self, cat, frame, scale=3.0, ec='cyan', cc='red'):
        """
        Draw hugs-pipe measurements.

        Parameters
        ----------
        cat : pandas.DataFrame or astropy.table.Table 
            Source catalog from hugs-pipe.
        frame : int
            ds9 frame (must exist)
        scale : float, optional
            Factor to multiply axes by. 
        ec : str, optional
            Ellipse color. 
        cc : str, optional
            Centroid color.
        """

        cat = check_astropy_to_pandas(cat)
        disp = self.frames[frame].disp
        bbox = self.frames[frame].bbox

        with disp.Buffering():
            for _, source in cat.iterrows():
                x, y = source['x_hsc'], source['y_hsc']
                point = lsst.afw.geom.Point2I(int(x), int(y))
                if bbox.contains(point):
                    a = source['semimajor_axis_sigma']
                    b = source['semiminor_axis_sigma']
                    theta = source['orientation']
                    ell = afwGeom.ellipses.Axes(scale*a, scale*b, theta)
                    disp.dot(ell, x, y, ctype=ec)
                    disp.dot('+', x, y, ctype=cc)

    def mpl_display_patch(self, hsc_bands='IRG', show=False, 
                          subplots_kw={}, imshow_kw={}, subplots=None,
                          rgb_kw={'Q':8, 'dataRange':1.25}):
        """
        Display RGB color image of patch using matplotlib.

        Parameters
        ----------

        Returns
        -------
        """

        images = {'R': None, 'G': None, 'B': None}
        for rgb_label, band in zip('RGB', hsc_bands): 
            images[rgb_label] = self.rgb_images[band]

        img = afwRgb.makeRGB(images['R'], 
                             images['G'],
                             images['B'],
                             **rgb_kw)

        if subplots is None:
            fig, ax = plt.subplots(
                subplot_kw={'xticks':[], 'yticks':[]}, **subplots_kw)
        else:
            fig, ax = subplots

        ax.imshow(img, origin='lower', **imshow_kw)

        if show:
            try: import RaiseWindow
            except: pass
            plt.show()

        return lsst.pipe.base.Struct(fig=fig, ax=ax)

    def mpl_display_cutout(self, coord_hsc, size=200, hsc_bands='IRG', 
                           show=False, subplots_kw={}, imshow_kw={}, 
                           subplots=None, rgb_kw={'Q':8, 'dataRange':1.25}):
        """
        Display RGB color image of cutout using matplotlib.

        Parameters 
        ----------

        Returns
        -------
        """

        cutout = {'R': None, 'G': None, 'B': None}
        for rgb_label, band in zip('RGB', hsc_bands): 
            rgb = self.rgb_images[band]
            mi = imtools.get_cutout(coord_hsc, size, exp=rgb)
            cutout[rgb_label] = mi

        img = afwRgb.makeRGB(cutout['R'], cutout['G'], cutout['B'], **rgb_kw)

        if subplots is None:
            fig, ax = plt.subplots(
                subplot_kw={'xticks':[], 'yticks':[]}, **subplots_kw)
        else:
            fig, ax = subplots

        ax.imshow(img, origin='lower', **imshow_kw)

        if show:
            try: import RaiseWindow
            except: pass
            plt.show()

        return lsst.pipe.base.Struct(fig=fig, ax=ax)

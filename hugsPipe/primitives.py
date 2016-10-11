from __future__ import division, print_function

import numpy as np
import lsst.afw.image
import lsst.afw.detection as afwDetect
from . import utils

__all__ = ['associate', 'image_threshold', 'postage_stamps']


def associate(mask, fpset, r_in=5, r_out=15, max_on_bit=20., 
              min_pix=1, plane_name='THRESH_HIGH', dilate=None):
    """
    Associate footprints in fpset with footprints in mask plane 
    'plane_name'. A footprint is associated with an object if 
    an on bit falls within an annulus centered with respect to 
    all its peaks.
    

    Parameters
    ----------
    mask : lsst.afw.image.imageLib.MaskU
        Mask object with plane named 'plane_name'.
    fpset : lsst.afw.detection.detectionLib.FootprintSet
        Set of footprints to associate with objects in mask.
    r_in : float, optional
        Inner radius in pixels of the association annulus.
    r_out : float, optional
        Outer radius in pixels of the association annulus.
    max_on_bit : int, optional
        Maximum number of on bits to consider as associated.
    min_pix : int, optional
        Footprints with less pixels will be masked. Set to 1 
        to ignore this option.
    plane_name : string, optional
        Name of the bit plane in mask to associate footprints with. 
    dilate : int, optional
        If not None, dilate association image with array of shape 
        (dilate, dilate). 

    Returns
    -------
    seg_assoc : 2D ndarray
        Segmentation image with non-zero values for all footprints 
        that are associated with an object in the mask. 

    Notes
    -----
    seg_assoc also includes footprints near objects with the 
    'BRIGHT_OBJECT' bit set.
    """
    x0, y0 = mask.getXY0()
    # False --> use footprint ids
    seg = fpset.insertIntoImage(False).getArray().copy()
    shape = seg.shape
    seg_assoc = np.zeros(shape, dtype=int)
    for foot in fpset.getFootprints():
        peaks = np.array([[p.getCentroid()[0]-x0, 
                           p.getCentroid()[1]-y0] for p in foot.getPeaks()])
        xc, yc = peaks.mean(axis=0)
        rows, cols = utils.annuli(yc, xc, r_in, r_out, shape=shape)
        ann_pix = mask.getArray()[rows, cols]
        on_bits = (ann_pix & mask.getPlaneBitMask(plane_name))!=0
        on_bits |= (ann_pix & mask.getPlaneBitMask('BRIGHT_OBJECT'))!=0
        if np.sum(on_bits) > max_on_bit:
            seg_assoc[seg==foot.getId()] = 1
        elif foot.getNpix() < min_pix:
            seg_assoc[seg==foot.getId()] = 1
    if dilate is not None:
        from scipy import ndimage
        dilator = np.ones((dilate, dilate))
        seg_assoc = ndimage.binary_dilation(seg_assoc, dilator)
        seg_assoc = seg_assoc.astype(int)
    return seg_assoc


def image_threshold(masked_image, thresh, thresh_type='stdev', npix=1, 
                    rgrow=None, isogrow=False, plane_name='', mask=None,
                    clear_mask=True):
    """
    Image thresholding. A bit mask will be set with name 'plane_name'.

    Parameters
    ----------
    masked_image : lsst.afw.image.imageLib.MaskedImageF
        A masked image object.
    thresh : float
        Threshold value.
    thresh_type : string, optional
        Threshold type: stdev, pixel_stdev, bitmask, value,
        or variace.
    npix : int, optional
        Minimum number of touching pixels in an object.
    rgrow : int, optional
        Number of pixels to grow footprints.
    isogrow : bool, optional
        If True, use (expensive) isotropic grow. 
    plane_name : string, optional
        Name of bit plane.
    mask : lsst.afw.image.imageLib.MaskU, optional
        Mask to set if not same as in masked_imaged
    clear_mask : bool, optional
        If True, clear the bit plane before thresholding

    Returns
    -------
    fp : lsst.afw.detection.detectionLib.FootprintSet
        Footprints assoicated with detected objects.
    """
    mask = masked_image.getMask() if mask is None else mask
    thresh_type = getattr(afwDetect.Threshold, thresh_type.upper())
    thresh = afwDetect.Threshold(thresh, thresh_type)
    fp = afwDetect.FootprintSet(masked_image, thresh, '', npix)
    if rgrow is not None:
        fp = afwDetect.FootprintSet(fp, rgrow, isogrow)
    if plane_name:
        mask.addMaskPlane(plane_name)
        if clear_mask:
            mask.clearMaskPlane(mask.getMaskPlane(plane_name))
        fp.setMask(mask, plane_name)
    return fp


def postage_stamps(exposure, fpset, grow=10):
    """
    Cutout postage stamps of given footprints. 

    Parameters
    ----------

    Returns
    -------
    """
    stamps = []
    for foot in fpset.getFootprints:
        bbox = foot.getBBox()
        if grow:
            bbox.grow(grow)
        exp_cutout = exposure.Factory(exposure, bbox, lsst.afw.image.PARENT) 
        stamps.append(exp_cutout)

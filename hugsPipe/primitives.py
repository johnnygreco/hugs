from __future__ import division, print_function

import numpy as np


__all__ = ['get_psf_sigma', 'associate', 'image_threshold']


def _annulus(row_c, col_c, r_in, r_out, shape):
    """
    Find the indices within an annulus embedded in 
    an array of the specified shape. 

    Parameters
    ----------
    row_c, col_c : float
        Ceentral coordinates of the annulus.
    r_in : float
        Inner radius of annulus. If r_in=0, the 
        annulus becomes a circle. 
    r_out : float
        Outer radius of annulus. 
    shape : tuple
        Shape of the host array.

    Returns
    -------
    row_idx, col_idx : ndarray of int
        Indices of annulus within host array.
    """
    # find bounding box of annulus
    center = np.array([row_c, col_c])
    ul = np.ceil(center - r_out).astype(int)
    lr = np.floor(center + r_out).astype(int)
    ul = np.maximum(ul, np.array([0, 0]))
    lr = np.minimum(lr, np.array(shape[:2]) - 1)
    bb_shape = lr - ul + 1

    # generate the coords within annulus in bbox
    ctr_shift = center - ul
    bb_row, bb_col = np.ogrid[0:float(bb_shape[0]), 0:float(bb_shape[1])]
    radius = np.sqrt((bb_row-ctr_shift[0])**2 + (bb_col-ctr_shift[1])**2)
    row_idx, col_idx = np.nonzero((radius >= r_in) & (radius < r_out))

    # shift back to original coords
    row_idx.flags.writeable = True
    col_idx.flags.writeable = True
    row_idx += ul[0]
    col_idx += ul[1]
    return row_idx, col_idx
    

def get_psf_sigma(exposure):
    """
    Get sigma of point-spread function.

    Parameters
    ----------
    exposure : lsst.afw.image.imageLib.ExposureF
        Exposure object. 

    Returns
    -------
    sigma : float
        Approximate standard deviation of the PSF
        of the image in this exposure.
    """
    psf = exposure.getPsf()
    sigma = psf.computeShape().getDeterminantRadius()
    return sigma


def associate(mask, fpset, r_in=5, r_out=15, max_on_bit=20., 
              plane_name='THRESH_HIGH'):
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
    plane_name : string, optional
        Name of the bit plane in mask to associate footprints with. 

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
    seg_assoc = fpset.insertIntoImage(False).getArray().copy() 
    for foot in fpset.getFootprints():
        peaks = np.array([[p.getCentroid()[0]-x0, 
                           p.getCentroid()[1]-y0] for p in foot.getPeaks()])
        xc, yc = peaks.mean(axis=0)
        rows, cols = _annulus(yc, xc, r_in, r_out, shape=mask.getArray().shape)
        ann_pix = mask.getArray()[rows, cols]
        on_bits = (ann_pix & mask.getPlaneBitMask(plane_name))!=0
        on_bits |= (ann_pix & mask.getPlaneBitMask('BRIGHT_OBJECT'))!=0
        if np.sum(on_bits)<max_on_bit:
            seg_assoc[seg_assoc==foot.getId()] = 0
    return seg_assoc


def image_threshold(masked_image, thresh, thresh_type='stdev', npix=1, 
                    rgrow=None, isogrow=False, plane_name='', mask=None,
                    clear_mask=True):
    """
    Image thresholding. As bit mask will be set with name 'plane_name'.

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
    import lsst.afw.detection as afwDetect

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


def viz(image, transparency=75, frame=1, 
        colors=[('THRESH_HIGH', 'magenta'), ('THRESH_LOW', 'yellow')]):
    import lsst.afw.display as afwDisplay
    disp = afwDisplay.Display(frame)
    for name, color in colors:
        disp.setMaskPlaneColor(name, color) 
    disp.setMaskTransparency(transparency)
    disp.mtv(image)
    return disp

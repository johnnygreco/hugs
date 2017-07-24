from __future__ import division, print_function


from lsst.afw.image import Mask
import lsst.afw.math as afwMath
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase

__all__ = ['SigmaClippedConfig', 'SigmaClippedTask', 'get_clipped_sig_task']


class SigmaClippedConfig(pexConfig.Config):
    """
    Configuration for sigma clipped statistics task.
    """
    badMaskPlanes = pexConfig.ListField(
        doc = 'Mask planes that, if set, indicate the associated pixel '
              'not be included when the calculating statistics.',
        dtype = str,
        default = ('EDGE', 'BRIGHT_OBJECT'),
    )
    numSigmaClip = pexConfig.Field(
        doc = 'Number of sigmas at which to clip data',
        dtype = float,
        default = 3.0,
    )
    numIter = pexConfig.Field(
        doc = 'Number of iterations of sigma clipping',
        dtype = int,
        default = 2,
    )


class SigmaClippedTask(pipeBase.Task):
    """
    Task for performing sigma clipped statistics on image.
    """
    ConfigClass = SigmaClippedConfig
    _DefaultName = 'SigmaClippedTask'

    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        self._badPixelMask = MaskU.getPlaneBitMask(self.config.badMaskPlanes)
        self._statsControl = afwMath.StatisticsControl()
        self._statsControl.setNumSigmaClip(self.config.numSigmaClip)
        self._statsControl.setNumIter(self.config.numIter)
        self._statsControl.setAndMask(self._badPixelMask)

    def run(self, maskedImage):
        statObj = afwMath.makeStatistics(
            maskedImage, 
            afwMath.MEANCLIP | afwMath.STDEVCLIP | afwMath.ERRORS,
            self._statsControl
        )
        mean, meanErr = statObj.getResult(afwMath.MEANCLIP)
        stdev, stdevErr = statObj.getResult(afwMath.STDEVCLIP)
        return pipeBase.Struct(
            mean = mean,
            meanErr = meanErr,
            stdev = stdev,
            stdevErr = stdevErr
        )


def get_clipped_sig_task(nsigma_clip=3.0, niter=10, 
                         bad_mask_planes=('EDGE', 'BRIGHT_OBJECT')):
    """
    Helper function to get sigma clipped statistics task.

    Parameters 
    ----------
    nsigma_clip : float, optional
        Number of sigmas at which to clip data.
    niter : int, optional
        Number of iterations of sigma clipping.
    bad_mask_planes : tuple of strings
        Mask planes to exclude from the stats calculation. 

    Returns
    -------
    task : hugs.stats.SigmaClippedTask
        Clipped sigma stats task.

    """
    config = SigmaClippedConfig(numSigmaClip=nsigma_clip, numIter=niter,
                                badMaskPlanes=bad_mask_planes)
    task = SigmaClippedTask(config)
    return task

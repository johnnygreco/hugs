
import os
import pytest
import lsst.afw.image
from .. import primitives as prim

def test_deblend_stamps():
    dataDIR = os.environ.get('TEST_DATA_DIR')
    fn = os.path.join(dataDIR, 'detexp-HSC-I-9348-7,6.fits')
    exposure = lsst.afw.image.ExposureF(fn)
    table = prim.deblend_stamps(exposure)
    return table

if __name__=='__main__':
    table = test_deblend_stamps()

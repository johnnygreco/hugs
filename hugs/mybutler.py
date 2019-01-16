import os
import lsst.afw.image

class PersonalButler(object):
    """
    My personal butler for testing the software offline. 
    """

    _DEFAULT_DATA_PATH='/Users/jgreco/local-io/hsc-test-data'

    def __init__(self, data_path=None): 

        self.data_path = data_path if data_path else self._DEFAULT_DATA_PATH

    def get(self, coadd_label, data_id, **kwargs):

        fn = 'calexp-{}-{}-{}.fits'.format(
            data_id['filter'], data_id['tract'], data_id['patch'])
        fn = os.path.join(self.data_path, fn)

        if not os.path.isfile(fn):
            raise Exception(fn + 'does not exist!')

        if 'filename' in coadd_label:
            return [fn]
        else:
            return lsst.afw.image.ExposureF(fn)


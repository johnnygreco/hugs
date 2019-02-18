from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import os.path as op
import subprocess

from ._settings import PARAMETERS      # parameter file 
from ._settings import SEX_CONFIG_DIR  # config directory

__all__ = ['Wrapper']

class Wrapper(object):
    """
    A python wrapper for SExtractor.

    Parameters
    ----------
    config : dict, optional
        sextractor config parameters 
    io_dir : str, optional
        sextractor input/output directory

    Notes
    -----
        This wrapper was written for SExtractor version 2.19.5.
    """
    
    def __init__(self, config={}, io_dir=''):

        self._config = config
        self._config['PARAMETERS_NAME'] = PARAMETERS
        self._io_dir = io_dir

    @property
    def config(self):
        """
        Return the current config.
        """
        return self._config

    def get_io_dir(self, join=None): 
        """
        Return the sextractor io directory.

        Parameters
        ----------
        join : string, optional
            File or directory to join to the path.
        """
        return op.join(self._io_dir, join) if join else self._io_dir

    def print_config(self):
        """
        Print the current SExtractor configuration
        """
        print('\ncurrent config', '\n--------------')
        for k, v in self.config.items():
            print('{:<16} {:<16}'.format(k, v))

    def set_check_images(self, which='a', prefix=''):
        """
        Set sextractor's CHECKIMAGE config params. 

        Parameters
        ----------
        which : string or list of strings, optional
            The type(s) of check images. Accepts the 
            full type name or the following short 
            abbreviations: 
            -----------------------
            b     |  BACKGROUND 
            -b    |  -BACKGROUND
            f     |  FILTERED
            s     |  SEGMENTATION 
            o     |  OBJECTS
            -o    |  -OBJECTS 
            brms  |  BACKGROUND_RMS
            mb    |  MINIBACKGROUND
            mbrms |  MINIBACK_RMS
            -----------------------
            If 'all' is given, then all check-images 
            will be saved.
        prefix : string, optional
            Prefix to add to the output image names, which 
            by default are image_type.fits
        """
        imgtypes = {'b':'BACKGROUND',
                    '-b':'-BACKGROUND',
                    'f':'FILTERED',
                    'a':'APERTURES',
                    's':'SEGMENTATION',
                    'o':'OBJECTS',
                    '-o':'-OBJECTS',
                    'brms':'BACKGROUND_RMS',
                    'mb':'MINIBACKGROUND',
                    'mbrms':'MINIBACK_RMS'}
        which = list(imgtypes.keys()) if which=='all' else which
        which = which if type(which) is list else [which]
        for val in list(imgtypes.values()):
            imgtypes.update({val:val})
        fn = [prefix+imgtypes[w].replace('-','bcksub-')+'.fits' for w in which]
        for i in range(len(fn)):
            fn[i] = self.get_io_dir(fn[i])
        checkimg_names = ','.join(fn)
        checkimg_type = ','.join([imgtypes[w] for w in which])
        self._config['CHECKIMAGE_NAME'] = checkimg_names
        self._config['CHECKIMAGE_TYPE'] = checkimg_type

    def write_config_file(self, fn='sex.config'):
        """
        Write the current configuration file. This can be used as input 
        for sextractor, but it is intended to save as a log.

        Paramerers
        ----------
        fn : string, optional
            Config file name. 
        """
        fn = op.join(self._io_dir, fn) 
        out = open(fn, 'w')
        for key, val in self.config.items():
            print('{:<16} {:<16}'.format(k, v))
        out.close()

    def run(self, img_fn, cat_fn=None):
        """
        Run sextractor.

        Parameters
        ----------
        img_fn : string, optional
            Image file name. The images must be in sexin/{relpath} directory, which 
            was defined upon instantiation. 
        cat_fn : string, optional
            Catalog name. 
        """

        cat_fn = cat_fn if cat_fn else self.get_io_dir('sex.cat') 
        default_config_fn = op.join(SEX_CONFIG_DIR, 'hugs-default.sex')

        cmd = 'sex -c {} {}'.format(default_config_fn, img_fn)
        cmd += ' -CATALOG_NAME '+cat_fn

        for key, val in self.config.items():
            cmd += ' -'+key+' '+str(val)

        print('\nrunning', '\n-------\n'+cmd+'\n')
        subprocess.call(cmd, shell=True)

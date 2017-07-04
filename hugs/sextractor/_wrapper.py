from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import os.path as op
import subprocess

from ._settings import DEFAULT_CONFIG  # Default SExtractor configuration
from ._settings import DEFAULT_PARAMS  # Default parameters 

__all__ = ['Wrapper']

class Wrapper(object):
    """
    A python wrapper for SExtractor.

    Parameters
    ----------
    config : dict, optional
        Default config parameters that will remain fixed 
        within this SexWrapper instance.
    params : list, optional
        List of catalog parameters to save. Will be appended to 
        default params.
    sexin: string, optional
        The input sextractor directory. Input images must
        be within this directory. If 'default', the directory
        structure is assumed to be SEX_IO_DIR/sexin.
    sexout: string, optional
        The input sextractor directory. Output files will
        be written to this directory, which will be created
        if it doesn't exist. If 'default', the directory structure 
        is assumed to be SEX_IO_DIR/sexout.
    configdir: string, optional
        Directory containing sextractor's config files. If 
        'default', set to the SEX_CONFIG_DIR env variable.
    relpath : string, optional
        Relative path from within the sexin and sexout directories. 
        This param is a bit pointless if you provide custom 
        paths for sexin and sexout. 

    Notes
    -----
     i) DEFAULT_CONFIG is for configuration settings that 
        are different from sextractor's default but rarely ever change. 
        These parameters are overridden if given in the config dict. 
    ii) This wrapper was written for SExtractor version 2.19.5.
    """
    
    def __init__(self, config=DEFAULT_CONFIG, params=DEFAULT_PARAMS, 
                 io_dir=''):

        self._config = config
        self._params = params
        self._io_dir = io_dir

        if 'PARAMETERS_NAME' not in config.keys():
            self._config['PARAMETERS_NAME'] = op.join(io_dir, 'sex.param')

    @property
    def params(self):
        """
        Return the current catalog parameters.
        """
        return self._params
    
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

    def print_params(self):
        """
        Print the current catalog parameters.
        """
        print('\ncurrent parameters', '\n------------------')
        print('\n'.join(self.params))

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
        use_outdir : bool, optional
            If True, save file in the current outdir.
        """
        fn = op.join(self._io_dir, fn) 
        out = open(fn, 'w')
        for key, val in self.config.items():
            print('{:<16} {:<16}'.format(k, v))
        out.close()

    def write_param_file(self):
        """
        Write the current parameters file. By default, the 
        file is written to the current config directory.

        Paramerers
        ----------
        fn : string, optional
            Parameter file name. If 'from_config', 
            use the same name as in config. 
        in_config_dir : bool, optional
            If True, write the param file in the config 
            directory. Otherwise, the desired path must 
            be given in fn. 
        """
        fn = op.join(self._io_dir, self._config['PARAMETERS_NAME'])
        param_file = open(fn, 'w')
        print('\n'.join(self.params), file=param_file)
        param_file.close()

    def run(self, img_fn, cat_fn=None):
        """
        Run sextractor.

        Parameters
        ----------
        img_fn : string, optional
            Image file name. The images must be in sexin/{relpath} directory, which 
            was defined upon instantiation. 
        configfile : string, optional
            Sextractor configure file name. Must be in the current config 
            directory. All config params set in this instance will override 
            settings in this file. 
        cat : string, optional
            Catalog name. Will be written to sexout/{relpath} directory.

        Notes
        -----
         i) All the sextractor configuration files must be in the config directory. 
        ii) Input files are assumed to be in sexin/{relpath}, and output files will
            be saved to sexout/{relpath}, which wilself._io_dir,l be created if it doesn't exist. 
        """

        cat_fn = cat_fn if cat_fn else self.get_io_dir('sex.cat') 

        cmd = 'sex '+img_fn
        cmd += ' -CATALOG_NAME '+cat_fn

        for key, val in self.config.items():
            cmd += ' -'+key+' '+str(val)

        self.write_param_file()
        print('\nrunning', '\n-------\n'+cmd+'\n')
        subprocess.call(cmd, shell=True)

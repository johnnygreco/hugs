from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import os.path as op

from ._settings import DEFAULT_CONFIG  # Default SExtractor configuration
from ._settings import DEFAULT_PARAMS  # Default catalog parameters to save
from ._settings import SEX_CONFIG_DIR
from ._settings import SEX_IO_DIR

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
    
    def __init__(self, config={}, params=[], config_dir=SEX_CONFIG_DIR):

        #################################################
        # set default configuration for this instance
        #################################################

        self._default_config = DEFAULT_CONFIG.copy()
        for key, val in config.items():
            self._default_config[key.upper()] = val
        self.reset_config()

        #################################################
        # If not given in config, write the catalog 
        # parameter file. 
        #################################################

        self._default_params = DEFAULT_PARAMS
        for par in params:
            if par not in DEFAULT_PARAMS:
                self._default_params.append(par)
        self.reset_params()
        self.write_param_file()

    def reset_config(self):
        """
        Reset the sextractor config params to the defaults 
        given upon instantiation. 
        """
        self._config = self._default_config.copy()

    def reset_params(self):
        """
        Reset the sextractor catalog params to the defaults 
        given upon instantiation. 
        """
        self._params = self._default_params

    def reset(self):
        """
        Reset the sextractor catalog & config params to the defaults 
        given upon instantiation. 
        """
        self.reset_config()
        self.reset_params()

    def get_sexin(self, join=None): 
        """
        Return the sexin input directory.

        Parameters
        ----------
        join : string, optional
            File or directory to join to the path.
        """
        return self._sexin if join is None else op.join(self._sexin, join)

    def get_sexout(self, join=None): 
        """
        Return the sexout output directory.

        Parameters
        ----------
        join : string, optional
            File or directory to join to the path.
        """
        return self._sexout if join is None else op.join(self._sexout, join)

    def get_indir(self, join=None): 
        """
        Return the full input directory path.
        If relpath not given, indir=sexin.

        Parameters
        ----------
        join : string, optional
            File or directory to join to the path.
        """
        return self._indir if join is None else op.join(self._indir, join)

    def get_outdir(self, join=None): 
        """
        Return the full output directory path.
        If relpath not given, outdir=sexout.

        Parameters
        ----------
        join : string, optional
            File or directory to join to the path.
        """
        return self._outdir if join is None else op.join(self._outdir, join)

    def get_configdir(self, join=None): 
        """
        Return the current config directory. 

        Parameters
        ----------
        join : string, optional
            File or directory to join to the path.
        """
        return self._configdir if join is None else op.join(self._configdir, join)

    def print_config(self):
        """
        Print the current SExtractor configuration
        """
        print('\ncurrent config', '\n--------------')
        for k, v in self.get_config().items():
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
        checkimg_names = ','.join(fn)
        checkimg_type = ','.join([imgtypes[w] for w in which])
        self._config['CHECKIMAGE_NAME'] = checkimg_names
        self._config['CHECKIMAGE_TYPE'] = checkimg_type

    def write_config_file(self, fn='myconfig.sex', use_outdir=False):
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
        fn = op.join(self._outdir, fn) if use_outdir else fn
        out = open(fn, 'w')
        for key, val in self.get_config().items():
            print('{:<16} {:<16}'.format(k, v))
        out.close()

    def add_param(self, param):
        """
        Add catalog parameter(s) to save.

        Parameters
        ----------
        param : string or list of strings
        """
        if type(param) is not list:
            param = [param]
        self._params += param
        self.write_param_file()

    def remove_param(self, param):
        """
        Remove catalog parameter(s).

        Parameters
        ----------
        param : string or list of strings
        """
        if type(param) is list:
            for p in param:
                self._params.remove(p)
        else:
            self._params.remove(param)
        self.write_param_file()

    def write_param_file(self, fn='from_config', in_config_dir=True):
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
        if fn=='from_config':
            fn = self._config['PARAMETERS_NAME']
        if in_config_dir:
            fn = op.join(self._configdir, fn)
        param_file = open(fn, 'w')
        print('\n'.join(self._params), file=param_file)
        param_file.close()

    def get_params(self):
        """
        Return the current catalog parameters.
        """
        return self._params

    def print_params(self):
        """
        Print the current catalog parameters.
        """
        print('\ncurrent parameters', '\n------------------')
        print('\n'.join(self.get_params()))

    def run(self, imagefile='img.fits', configfile='default.sex', cat='sex.cat'):
        """
        Run sextractor.

        Parameters
        ----------
        imagefile : string, optional
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
            be saved to sexout/{relpath}, which will be created if it doesn't exist. 
        """
        import subprocess

        # get the run directory 
        rundir = os.getcwd()

        # run from config directory
        os.chdir(self._configdir)

        if ',' in imagefile:
            imfile = imagefile.split(',')
            assert len(imfile)==2
            imfile = op.join(self._indir, imfile[0])+','+\
                     op.join(self._indir, imfile[1])
        else:
            imfile = op.join(self._indir, imagefile)
        catfile = op.join(self._outdir, cat)
        cmd = 'sex -c '+configfile+' '+imfile
        cmd += ' -CATALOG_NAME '+catfile

        for key, val in self._config.items():
            if 'IMAGE' in key:
                if 'CHECK' not in key:
                    val = op.join(self._indir, val)
                elif key == 'CHECKIMAGE_NAME':
                    withpath = [op.join(self._outdir,v) for v in val.split(',')]
                    val = ','.join(withpath)
            if key=='ASSOC_NAME':
                val = op.join(self._outdir, val)
            cmd += ' -'+key+' '+str(val)

        print('\nrunning', '\n-------\n'+cmd+'\n')
        subprocess.call(cmd, shell=True)

        # change back to original run directory
        os.chdir(rundir)

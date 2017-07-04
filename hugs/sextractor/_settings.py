from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from ..utils import project_dir

__all__ = ['SEX_CONFIG_DIR', 'DEFAULT_CONFIG', 'DEFAULT_PARAMS']

# My default location for sextractor's config files
SEX_CONFIG_DIR = os.path.join(project_dir, 'hugs/sextractor/config')

# Default sextractor configuration 
default_config_file = os.path.join(SEX_CONFIG_DIR, 'hugs-sextractor.config')
if os.path.isfile(default_config_file):
    DEFAULT_CONFIG = {}
    with open(default_config_file) as f:
        l = f.read().split()
    for i in range(0, len(l), 2):
        if (l[i]=='PARAMETERS_NAME') or (l[i]=='FILTER_NAME'):
            DEFAULT_CONFIG[l[i]] = os.path.join(SEX_CONFIG_DIR, l[i+1])
        else:
            DEFAULT_CONFIG.update({l[i]:l[i+1]})

# Default catalog parameters. 
default_param_file = os.path.join(SEX_CONFIG_DIR, 'hugs-sextractor.param')
if os.path.isfile(default_param_file):
    with open(default_param_file) as f:
        DEFAULT_PARAMS = f.read().split()

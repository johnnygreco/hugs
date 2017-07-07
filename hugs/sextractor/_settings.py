from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
from ..utils import project_dir

__all__ = ['SEX_CONFIG_DIR', 'PARAMETERS']

# My default location for sextractor's config files
SEX_CONFIG_DIR = os.path.join(project_dir, 'hugs/sextractor/config')
PARAMETERS = os.path.join(
    project_dir, 'hugs/sextractor/config/hugs-sextractor.param')

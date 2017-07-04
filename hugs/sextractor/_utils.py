from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os

__all__ = ['read_cat', 'cat_to_ds9reg']

def read_cat(catfile):
    """
    Read SExtractor output catalog using astropy.

    Parameters
    ----------
    catfile : string
      SExtractor output catalog.

    Returns
    -------
    cat : astropy table
      Output sextractor catalog.
    """
    from astropy.table import Table
    if os.stat(catfile).st_size==0:
        cat = Table()
    else:
        from collections import OrderedDict
        cat = Table.read(catfile, format='ascii.sextractor')
        cat.meta = OrderedDict() # bug in astropy sextractor table
    return cat

def cat_to_ds9reg(cat, color='green', tag='all', winparams=False,
                  outfile='sex.reg', drawmode='ellipse', textparam=None):
    """
    Write a ds9 region file from SExtractor catalog. 

    Parameters
    ----------
    cat : astropy.table.Table
        Output cat from a sextractor run (output of read_cat).
    color : string, optional 
      Region color (cyan, blue, magenta, red, green, 
      yellow, white, or black)
    tag : string, optional
      ds9 tag for all the regions
    winparams : bool, optional
        If True, use sextractor's windowed parameters.
    outfile : string, optional
        Output reg file name. 
    drawmode : string, optional
        Draw an 'ellipse' or 'point' for every object
    textparam : string, optional
        If not None, write this sextractor output parameter 
        next to each object in the catalog. 

    Notes
    -----
     i) Adapted from https://github.com/nhmc/Barak.git. 
    ii) The sextractor output file must contain X_IMAGE and 
        Y_IMAGE for drawmode=point, and for drawmode=ellipse, 
        it must also contain A_IMAGE, B_IMAGE, and THETA_IMAGE.
        The corresponding 'WIN' parameters are acceptable with 
        winparams set to True. 
    """
    assert (drawmode=='ellipse') or (drawmode=='point')

    regions = ['global font="helvetica 10 normal" select=1 highlite=1 '
               'edit=0 move=1 delete=1 include=1 fixed=0 source']
    regions.append('image')

    fields = ['X_IMAGE', 'Y_IMAGE', 'A_IMAGE','B_IMAGE','THETA_IMAGE']
    if drawmode=='point':
        fields = fields[:2]
    if winparams:
        fields = [f.split('_')[0]+'WIN'+'_'+f.split('_')[1] for f in fields]
    if textparam is not None:
        textfmt = 'text={%s}'
        fields.append(textparam)
    else:
        textfmt = ''

    fmt = {'ellipse':'ellipse(%s %s %s %s %s) # '+textfmt+' color=%s tag={%s}',
           'point':'point(%s %s) # point=circle '+textfmt+' color=%s tag={%s}'}[drawmode]

    for row  in cat[fields]:
        vals = list(row)
        vals.extend([color, tag])
        regions.append(fmt % tuple(vals))

    print('writing to region file to', outfile)
    fh = open(outfile,'w')
    fh.write('\n'.join(regions))
    fh.close()

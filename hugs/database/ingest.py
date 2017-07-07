from __future__ import division, print_function


from sqlalchemy import exists
from sqlalchemy.sql import func, and_

from .tables import Run, Tract, Patch
from .tables import Source, AperturePhotometry, CircularAperture
from .connect import connect, Session

__all__ = ['HugsIngest', 'add_patch_catalog']

def _get_current_id(session, table_id):
    return session.query(func.max(table_id)).first()[0]

class HugsIngest(object):

    def __init__(self, session, run_name):

        self.run_name = run_name 
        self.session = session
        self.current_tract_id = None
        self.current_patch_id = None
        run_query = self.session.query(Run).filter(Run.name==run_name)
        num_rows = run_query.count()
        if num_rows==0:
            self.session.add(Run(name=run_name))
            self.session.commit()
            self.run_id =  _get_current_id(self.session, Run.id)
        elif num_rows==1:
            self.run_id = run_query.first().id
        else:
            print('Warning {} name {}... db needs attention'.format(
                num_rows, run_name))
        
    def add_tract(self, tract):
        tract_query = self.session.query(Tract).filter(
            and_(Tract.hsc_id==tract, Tract.run_id==self.run_id))
        num_rows = tract_query.count()
        if num_rows==0:
            self.session.add(Tract(hsc_id=tract, run_id=self.run_id))
            self.session.commit()
            self.current_tract_id =  _get_current_id(self.session, Tract.id)
        elif num_rows==1:
            self.current_tract_id = tract_query.first().id
        else:
            print('Warning {} rows with tract {}... db needs attention'.format(
                num_rows, tract))

    def add_patch(self, patch, patch_meta):
        patch_row = Patch(
            hsc_id=patch, 
            x0=patch_meta.x0, 
            y0=patch_meta.y0, 
            cleaned_frac=patch_meta.cleaned_frac, 
            bright_obj_frac=patch_meta.bright_obj_frac,
            good_data_frac=patch_meta.good_data_frac,
            tract_id=self.current_tract_id
        )
        self.session.add(patch_row) 
        self.session.commit()
        self.current_patch_id = _get_current_id(self.session, Patch.id)

    def add_catalog(self, catalog, aper_radii=[]):
        """
        """
        assert self.current_tract_id is not None
        assert self.current_patch_id is not None

        # ingest source catalog
        sources = []
        aper_phot = []
        circ_aper = []
        aper_phot_id = 0
        for row, obj in enumerate(catalog):
            source_id = row + 1
            src = Source(
                x_image=obj['X_IMAGE'], 
                y_image=obj['Y_IMAGE'], 
                ra=obj['ALPHA_J2000'], 
                dec=obj['DELTA_J2000'], 
                a_image=obj['A_IMAGE'], 
                b_image=obj['B_IMAGE'],
                theta_image=obj['THETA_IMAGE'], 
                ellipticity=obj['ELLIPTICITY'],
                kron_radius=obj['KRON_RADIUS'], 
                patch_id=self.current_patch_id
            )
            sources.append(src)
            aper_phot_id += 1
            for band in 'gri':
                phot = AperturePhotometry(
                    bandpass=band, 
                    mag_auto=obj['MAG_AUTO('+band+')'],
                    fwhm_image=obj['FWHM_IMAGE('+band+')'], 
                    flux_radius=obj['FLUX_RADIUS('+band+')'], 
                    source_id=source_id
                )
                aper_phot.append(phot)
                for num, rad in enumerate(aper_radii):
                    circ = CircularAperture(
                        mag=obj['MAG_APER_{}({})'.format(num, band)], 
                        radius=rad, 
                        aper_phot_id=aper_phot_id
                    )
                    circ_aper.append(circ)
        self.session.add_all(sources)
        self.session.add_all(aper_phot)
        self.session.add_all(circ_aper)
        self.session.commit()

    def add_all(self, tract, patch, patch_meta, catalog, aper_radii=[]):
        self.add_tract(tract)
        self.add_patch(patch, patch_meta)
        self.add_catalog(catalog, aper_radii)

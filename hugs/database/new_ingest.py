from __future__ import division, print_function


from sqlalchemy import exists
from sqlalchemy.sql import func, and_

from .new_tables import Run, Tract, Patch, Source
from .connect import connect, Session
from ..utils import pixscale, ext_coeff, get_dust_map
dustmap = get_dust_map()

__all__ = ['HugsIngest']


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
            self.run_id =  self._get_current_id(Run.id)
        elif num_rows==1:
            self.run_id = run_query.first().id
        else:
            print('Warning {} name {}... db needs attention'.format(
                num_rows, run_name))

    def _get_current_id(self, table_id):
        return self.session.query(func.max(table_id)).first()[0]
        
    def add_tract(self, tract):
        tract_query = self.session.query(Tract).filter(
            and_(Tract.hsc_id==tract, Tract.run_id==self.run_id))
        num_rows = tract_query.count()
        if num_rows==0:
            self.session.add(Tract(hsc_id=tract, run_id=self.run_id))
            self.session.commit()
            self.current_tract_id = self._get_current_id(Tract.id)
        elif num_rows==1:
            self.current_tract_id = tract_query.first().id
        else:
            print('Warning {} rows with tract {}... db needs attention'.format(
                num_rows, tract))

    def add_patch(self, patch, patch_meta):
        assert self.current_tract_id is not None
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
        self.current_patch_id = self._get_current_id(Patch.id)

    def add_catalog(self, catalog, num_apertures=5):
        """
        """
        assert self.current_patch_id is not None

        # ingest source catalog
        sources = []
        for i, obj in enumerate(catalog):
            src = dict(
                x=obj['x_img'], 
                y=obj['y_img'], 
                ra=obj['ALPHA_J2000'], 
                dec=obj['DELTA_J2000'], 
                a_image=obj['A_IMAGE'], 
                b_image=obj['B_IMAGE'],
                theta_image=obj['THETA_IMAGE'], 
                ellipticity=obj['ELLIPTICITY'],
                kron_radius=obj['KRON_RADIUS'], 
                petro_radius=obj['PETRO_RADIUS'], 
                flags=obj['FLAGS'],
                patch_id=self.current_patch_id
            )
            for b in 'gri':
                ebv = dustmap.ebv(obj['ALPHA_J2000'], obj['DELTA_J2000'])
                A_lam = ebv*getattr(ext_coeff, b)
                src.update({
                    'mag_auto_'+b : obj['MAG_AUTO('+b+')'],
                    'mag_auto_'+b+'_err' : obj['MAGERR_AUTO('+b+')'],
                    'mag_petro_'+b : obj['MAG_PETRO('+b+')'],
                    'mag_petro_'+b+'_err' : obj['MAGERR_PETRO('+b+')'],
                    'fwhm_'+b : obj['FWHM_IMAGE('+b+')']*pixscale, 
                    'flux_radius_'+b : obj['FLUX_RADIUS('+b+')']*pixscale, 
                    'ebv' : ebv,
                    'A_'+b : dustmap.ebv(obj['ALPHA_J2000'])
                })
                for num in range(num_apertures):
                    mag_ap = 'mag_ap{}_{}'.format(num, b)
                    mag_ap_err = 'mag_ap{}_{}_err'.format(num, b)
                    src.update({
                        mag_ap : obj['MAG_APER_{}({})'.format(num, b)], 
                        mag_ap_err : obj['MAGERR_APER_{}({})'.format(num, b)], 
                    })
            sources.append(src)
            if i % 10 == 0 :
                self.session.execute(Source.__table__.insert(), sources)
                self.session.commit()
                sources = []
        if len(sources)>0:
            self.session.execute(Source.__table__.insert(), sources)
            self.session.commit()

    def add_all(self, tract, patch, patch_meta, catalog, num_apertures=5): 
        self.add_tract(tract)
        self.add_patch(patch, patch_meta)
        self.add_catalog(catalog, num_apertures)

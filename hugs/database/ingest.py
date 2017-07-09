from __future__ import division, print_function


from sqlalchemy import exists
from sqlalchemy.sql import func, and_

from .tables import Run, Tract, Patch, Source
from .connect import connect, Session

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
            print('Warning {} rows in run name {}'.format(num_rows, run_name))

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
            print('Warning {} rows with tract {}'.format(num_rows, tract))

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

    def add_catalog(self, catalog):
        """
        """
        assert self.current_patch_id is not None
        catalog['patch_id'] = self.current_patch_id
        catalog.to_sql('source', self.session.bind, 
                       if_exists='append', index=False)

    def add_all(self, tract, patch, patch_meta, catalog): 
        self.add_tract(tract)
        self.add_patch(patch, patch_meta)
        self.add_catalog(catalog)

from __future__ import division, print_function

import numpy as np
from sqlalchemy import String, Integer, Float, Column
from sqlalchemy.schema import ForeignKey
from sqlalchemy.orm import relationship
from .new_connect import Base
from astropy.coordinates import SkyCoord
from ..utils import pixscale, ext_coeff, get_dust_map
dustmap = get_dust_map()

__all__ = ['Run', 'Tract', 'Patch', 'Source']


class Run(Base):
    __tablename__ = 'run'

    # Table columns
    id = Column(Integer, primary_key=True)
    name = Column(String, nullable=False)

    # Relationships
    Tracts = relationship('Tract', cascade='all, delete-orphan')


class Tract(Base):
    __tablename__ = 'tract'

    # Table columns
    id = Column(Integer, primary_key=True)
    hsc_id = Column(Integer, nullable=False)

    # Relationships
    run_id = Column(Integer, ForeignKey('run.id'), nullable=False)
    run = relationship('Run')
    patches = relationship('Patch', cascade='all, delete-orphan')


class Patch(Base):
    __tablename__ = 'patch'

    # Table columns
    id = Column(Integer, primary_key=True)
    hsc_id = Column(String, nullable=False) 

    x0 = Column(Float, nullable=False)
    y0 = Column(Float, nullable=False)

    good_data_frac = Column(Float, nullable=True)
    cleaned_frac = Column(Float, nullable=True)
    bright_obj_frac = Column(Float, nullable=True)

    # Relationships
    tract_id = Column(Integer, ForeignKey('tract.id'), nullable=False)
    tract = relationship('Tract')
    sources = relationship('Source', cascade='all, delete-orphan')


class Source(Base):
    __tablename__ = 'source'
    
    # Table columns
    id = Column(Integer, primary_key=True)
    x = Column(Float, nullable=False)
    y = Column(Float, nullable=False)
    ra =  Column(Float, nullable=False)
    dec = Column(Float, nullable=False)

    a_image = Column(Float, nullable=True)  
    b_image = Column(Float, nullable=True)  
    theta_image = Column(Float, nullable=True)  
    ellipticity = Column(Float, nullable=True)  
    kron_radius = Column(Float, nullable=True)  
    petro_radius = Column(Float, nullable=True)  
    flags = Column(Integer, nullable=False)  

    mag_auto_g = Column(Float, nullable=True)
    mag_auto_r = Column(Float, nullable=True)
    mag_auto_i = Column(Float, nullable=True)
    mag_auto_g_err = Column(Float, nullable=True)
    mag_auto_r_err = Column(Float, nullable=True)
    mag_auto_i_err = Column(Float, nullable=True)

    mag_petro_g = Column(Float, nullable=True)
    mag_petro_r = Column(Float, nullable=True)
    mag_petro_i = Column(Float, nullable=True)
    mag_petro_g_err = Column(Float, nullable=True)
    mag_petro_r_err = Column(Float, nullable=True)
    mag_petro_i_err = Column(Float, nullable=True)

    mag_ap0_g = Column(Float, nullable=True)
    mag_ap1_g = Column(Float, nullable=True)
    mag_ap2_g = Column(Float, nullable=True)
    mag_ap3_g = Column(Float, nullable=True)
    mag_ap4_g = Column(Float, nullable=True)
    mag_ap0_g_err = Column(Float, nullable=True)
    mag_ap1_g_err = Column(Float, nullable=True)
    mag_ap2_g_err = Column(Float, nullable=True)
    mag_ap3_g_err = Column(Float, nullable=True)
    mag_ap4_g_err = Column(Float, nullable=True)

    mag_ap0_r = Column(Float, nullable=True)
    mag_ap1_r = Column(Float, nullable=True)
    mag_ap2_r = Column(Float, nullable=True)
    mag_ap3_r = Column(Float, nullable=True)
    mag_ap4_r = Column(Float, nullable=True)
    mag_ap0_r_err = Column(Float, nullable=True)
    mag_ap1_r_err = Column(Float, nullable=True)
    mag_ap2_r_err = Column(Float, nullable=True)
    mag_ap3_r_err = Column(Float, nullable=True)
    mag_ap4_r_err = Column(Float, nullable=True)

    mag_ap0_i = Column(Float, nullable=True)
    mag_ap1_i = Column(Float, nullable=True)
    mag_ap2_i = Column(Float, nullable=True)
    mag_ap3_i = Column(Float, nullable=True)
    mag_ap4_i = Column(Float, nullable=True)
    mag_ap0_i_err = Column(Float, nullable=True)
    mag_ap1_i_err = Column(Float, nullable=True)
    mag_ap2_i_err = Column(Float, nullable=True)
    mag_ap3_i_err = Column(Float, nullable=True)
    mag_ap4_i_err = Column(Float, nullable=True)

    fwhm_g = Column(Float, nullable=True)  
    fwhm_r = Column(Float, nullable=True)  
    fwhm_i = Column(Float, nullable=True)  
    flux_radius_g = Column(Float, nullable=True)
    flux_radius_r = Column(Float, nullable=True)
    flux_radius_i = Column(Float, nullable=True)

    ebv = Column(Float, nullable=True)
    A_g = Column(Float, nullable=True)
    A_r = Column(Float, nullable=True)
    A_i = Column(Float, nullable=True)

    # Relationships
    patch_id = Column(Integer, ForeignKey('patch.id'), nullable=False)
    patch = relationship('Patch')

    @property
    def x_hsc(self):
        return self.x + self.patch.x0

    @property
    def y_hsc(self):
        return self.y + self.patch.y0

    @property
    def skycoord(self):
        return SkyCoord(ra=self.ra, dec=self.dec, unit='deg')

    @property
    def hr_angle_string(self):
        return self.skycoord.to_string('hmsdms')

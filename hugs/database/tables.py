from __future__ import division, print_function

import numpy as np
from sqlalchemy import String, Integer, Float, Column
from sqlalchemy.schema import ForeignKey
from sqlalchemy.orm import relationship
from .connect import Base
from astropy.coordinates import SkyCoord

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
    small_frac = Column(Float, nullable=True)
    cleaned_frac = Column(Float, nullable=True)
    bright_obj_frac = Column(Float, nullable=True)

    # TODO
    # add number of detected sources column

    # Relationships
    tract_id = Column(Integer, ForeignKey('tract.id'), nullable=False)
    tract = relationship('Tract')
    sources = relationship('Source', cascade='all, delete-orphan')


class Source(Base):
    __tablename__ = 'source'

    # Table columns
    id = Column(Integer, primary_key=True)
    x_image = Column(Float, nullable=False)
    y_image = Column(Float, nullable=False)
    x_hsc = Column(Float, nullable=False)
    y_hsc = Column(Float, nullable=False)
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
    magerr_auto_g = Column(Float, nullable=True)
    magerr_auto_r = Column(Float, nullable=True)
    magerr_auto_i = Column(Float, nullable=True)

    mag_petro_g = Column(Float, nullable=True)
    mag_petro_r = Column(Float, nullable=True)
    mag_petro_i = Column(Float, nullable=True)
    magerr_petro_g = Column(Float, nullable=True)
    magerr_petro_r = Column(Float, nullable=True)
    magerr_petro_i = Column(Float, nullable=True)

    mag_ap0_g = Column(Float, nullable=True)
    mag_ap1_g = Column(Float, nullable=True)
    mag_ap2_g = Column(Float, nullable=True)
    mag_ap3_g = Column(Float, nullable=True)
    mag_ap4_g = Column(Float, nullable=True)
    mag_ap5_g = Column(Float, nullable=True)
    mag_ap6_g = Column(Float, nullable=True)
    mag_ap7_g = Column(Float, nullable=True)
    mag_ap8_g = Column(Float, nullable=True)
    mag_ap9_g = Column(Float, nullable=True)
    magerr_ap0_g = Column(Float, nullable=True)
    magerr_ap1_g = Column(Float, nullable=True)
    magerr_ap2_g = Column(Float, nullable=True)
    magerr_ap3_g = Column(Float, nullable=True)
    magerr_ap4_g = Column(Float, nullable=True)
    magerr_ap5_g = Column(Float, nullable=True)
    magerr_ap6_g = Column(Float, nullable=True)
    magerr_ap7_g = Column(Float, nullable=True)
    magerr_ap8_g = Column(Float, nullable=True)
    magerr_ap9_g = Column(Float, nullable=True)

    mag_ap0_r = Column(Float, nullable=True)
    mag_ap1_r = Column(Float, nullable=True)
    mag_ap2_r = Column(Float, nullable=True)
    mag_ap3_r = Column(Float, nullable=True)
    mag_ap4_r = Column(Float, nullable=True)
    mag_ap5_r = Column(Float, nullable=True)
    mag_ap6_r = Column(Float, nullable=True)
    mag_ap7_r = Column(Float, nullable=True)
    mag_ap8_r = Column(Float, nullable=True)
    mag_ap9_r = Column(Float, nullable=True)
    magerr_ap0_r = Column(Float, nullable=True)
    magerr_ap1_r = Column(Float, nullable=True)
    magerr_ap2_r = Column(Float, nullable=True)
    magerr_ap3_r = Column(Float, nullable=True)
    magerr_ap4_r = Column(Float, nullable=True)
    magerr_ap5_r = Column(Float, nullable=True)
    magerr_ap6_r = Column(Float, nullable=True)
    magerr_ap7_r = Column(Float, nullable=True)
    magerr_ap8_r = Column(Float, nullable=True)
    magerr_ap9_r = Column(Float, nullable=True)

    mag_ap0_i = Column(Float, nullable=True)
    mag_ap1_i = Column(Float, nullable=True)
    mag_ap2_i = Column(Float, nullable=True)
    mag_ap3_i = Column(Float, nullable=True)
    mag_ap4_i = Column(Float, nullable=True)
    mag_ap5_i = Column(Float, nullable=True)
    mag_ap6_i = Column(Float, nullable=True)
    mag_ap7_i = Column(Float, nullable=True)
    mag_ap8_i = Column(Float, nullable=True)
    mag_ap9_i = Column(Float, nullable=True)
    magerr_ap0_i = Column(Float, nullable=True)
    magerr_ap1_i = Column(Float, nullable=True)
    magerr_ap2_i = Column(Float, nullable=True)
    magerr_ap3_i = Column(Float, nullable=True)
    magerr_ap4_i = Column(Float, nullable=True)
    magerr_ap5_i = Column(Float, nullable=True)
    magerr_ap6_i = Column(Float, nullable=True)
    magerr_ap7_i = Column(Float, nullable=True)
    magerr_ap8_i = Column(Float, nullable=True)
    magerr_ap9_i = Column(Float, nullable=True)

    fwhm_g = Column(Float, nullable=True)  
    fwhm_r = Column(Float, nullable=True)  
    fwhm_i = Column(Float, nullable=True)  

    flux_radius_10_g = Column(Float, nullable=True)
    flux_radius_20_g = Column(Float, nullable=True)
    flux_radius_30_g = Column(Float, nullable=True)
    flux_radius_40_g = Column(Float, nullable=True)
    flux_radius_50_g = Column(Float, nullable=True)
    flux_radius_60_g = Column(Float, nullable=True)
    flux_radius_65_g = Column(Float, nullable=True)
    flux_radius_70_g = Column(Float, nullable=True)
    flux_radius_80_g = Column(Float, nullable=True)
    flux_radius_90_g = Column(Float, nullable=True)

    flux_radius_10_r = Column(Float, nullable=True)
    flux_radius_20_r = Column(Float, nullable=True)
    flux_radius_30_r = Column(Float, nullable=True)
    flux_radius_40_r = Column(Float, nullable=True)
    flux_radius_50_r = Column(Float, nullable=True)
    flux_radius_60_r = Column(Float, nullable=True)
    flux_radius_65_r = Column(Float, nullable=True)
    flux_radius_70_r = Column(Float, nullable=True)
    flux_radius_80_r = Column(Float, nullable=True)
    flux_radius_90_r = Column(Float, nullable=True)

    flux_radius_10_i = Column(Float, nullable=True)
    flux_radius_20_i = Column(Float, nullable=True)
    flux_radius_30_i = Column(Float, nullable=True)
    flux_radius_40_i = Column(Float, nullable=True)
    flux_radius_50_i = Column(Float, nullable=True)
    flux_radius_60_i = Column(Float, nullable=True)
    flux_radius_65_i = Column(Float, nullable=True)
    flux_radius_70_i = Column(Float, nullable=True)
    flux_radius_80_i = Column(Float, nullable=True)
    flux_radius_90_i = Column(Float, nullable=True)

    gini_full = Column(Float, nullable=True)
    gini_1 = Column(Float, nullable=True)
    gini_1p5 = Column(Float, nullable=True)
    gini_2 = Column(Float, nullable=True)
    gini_1p5_circ = Column(Float, nullable=True)
    gini_2_circ = Column(Float, nullable=True)
    
    acorr_peak = Column(Float, nullable=True)
    acorr_bkgd = Column(Float, nullable=True)
    acorr_ratio = Column(Float, nullable=True)

    ebv = Column(Float, nullable=True)
    A_g = Column(Float, nullable=True)
    A_r = Column(Float, nullable=True)
    A_i = Column(Float, nullable=True)

    # Relationships
    patch_id = Column(Integer, ForeignKey('patch.id'), nullable=False)
    patch = relationship('Patch')

    @property
    def skycoord(self):
        return SkyCoord(ra=self.ra, dec=self.dec, unit='deg')

    @property
    def hr_angle_string(self):
        return self.skycoord.to_string('hmsdms')

    
class Synth(Base):
    __tablename__ = 'synth'
    
    # Table columns
    id = Column(Integer, primary_key=True)
    synth_id = Column(Integer, nullable=False)
    mask_bright_object = Column(Integer, nullable=True)
    mask_cleaned = Column(Integer, nullable=True)
    mask_small = Column(Integer, nullable=True)
    mask_no_data = Column(Integer, nullable=True)
    mask_sat = Column(Integer, nullable=True)
    mask_suspect = Column(Integer, nullable=True)

    # Relationships
    patch_id = Column(Integer, ForeignKey('patch.id'), nullable=False)
    patch = relationship('Patch')

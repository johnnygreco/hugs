from __future__ import division, print_function

import numpy as np
from sqlalchemy import String, Integer, Float, Column
from sqlalchemy.schema import ForeignKey
from sqlalchemy.orm import relationship
from .connect import Base
from astropy.coordinates import SkyCoord
from ..utils import pixscale, ext_coeff, get_dust_map
dustmap = get_dust_map()

__all__ = ['Run', 'Tract', 'Patch', 'Source', 
           'AperturePhotometry', 'CircularAperture']


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

    # Relationships
    patch_id = Column(Integer, ForeignKey('patch.id'), nullable=False)
    patch = relationship('Patch')
    aper_phot = relationship('AperturePhotometry', 
                             cascade='all, delete-orphan')

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
    def ebv(self):
        return dustmap.ebv(self.ra, self.dec)

    @property
    def hr_angle_string(self):
        return self.skycoord.to_string('hmsdms')

    def get_band_phot(self, band):
        for phot in self.aper_phot:
            if phot.bandpass==band:
                return phot


class AperturePhotometry(Base):

    __tablename__ = 'aper_phot'

    # Table columns
    id  = Column(Integer, primary_key=True)
    bandpass = Column(String, nullable=False)
    mag_auto = Column(Float, nullable=True)
    mag_auto_err = Column(Float, nullable=True)
    mag_petro = Column(Float, nullable=True)
    mag_petro_err = Column(Float, nullable=True)
    fwhm_image = Column(Float, nullable=True)  
    flux_radius = Column(Float, nullable=True)

    # Relationships
    source_id = Column(Integer, ForeignKey('source.id'), nullable=False)
    source = relationship('Source')
    circ_aper = relationship('CircularAperture', 
                             cascade='all, delete-orphan')

    @property
    def flux_radius_arcsec(self):
        return self.flux_radius*pixscale

    @property
    def fwhm_arcsec(self):
        return self.fwhm_image*pixscale

    @property
    def A_lam(self):
        return getattr(ext_coeff, self.bandpass)*self.source.ebv


class CircularAperture(Base):

    __tablename__ = 'circ_aper'

    # Table columns
    id  = Column(Integer, primary_key=True)
    mag = Column(Float, nullable=True)
    mag_err = Column(Float, nullable=True)
    radius = Column(Float, nullable=False)

    # Relationships
    aper_phot_id = Column(Integer, ForeignKey('aper_phot.id'), nullable=False)
    aper_phot = relationship('AperturePhotometry')

    @property
    def area(self):
        return np.pi*self.radius**2

    @property
    def mu(self):
        return self.mag + 2.5*np.log10(self.area) 

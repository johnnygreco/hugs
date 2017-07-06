from __future__ import division, print_function

import numpy as np
from sqlalchemy import String, Integer, Float, Column
from sqlalchemy.schema import ForeignKey
from sqlalchemy.orm import relationship
from .connect import Base
from astropy.coordinates import SkyCoord

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
    x_image = Column(Float, nullable=False)
    y_image = Column(Float, nullable=False)
    ra =  Column(Float, nullable=False)
    dec = Column(Float, nullable=False)
    a_image = Column(Float, nullable=False)  
    b_image = Column(Float, nullable=False)  
    theta_image = Column(Float, nullable=False)  
    ellipticity = Column(Float, nullable=False)  
    kron_radius = Column(Float, nullable=False)  

    # Relationships
    patch_id = Column(Integer, ForeignKey('patch.id'), nullable=False)
    patch = relationship('Patch')
    aper_phot = relationship('AperturePhotometry', 
                             cascade='all, delete-orphan')

    @property
    def x_hsc(self):
        return self.x_image + self.patch.x0

    @property
    def y_hsc(self):
        return self.y_image + self.patch.y0

    @property
    def skycoord(self):
        return SkyCoord(ra=self.ra, dec=self.dec, unit='deg')

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
    mag_auto = Column(Float, nullable=False)
    fwhm_image = Column(Float, nullable=False)  

    # Relationships
    source_id = Column(Integer, ForeignKey('source.id'), nullable=False)
    source = relationship('Source')
    circ_aper = relationship('CircularAperture', 
                             cascade='all, delete-orphan')


class CircularAperture(Base):

    __tablename__ = 'circ_aper'

    # Table columns
    id  = Column(Integer, primary_key=True)
    mag = Column(Float, nullable=False)
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

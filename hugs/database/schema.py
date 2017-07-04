from __future__ import division, print_function

from sqlalchemy import Integer, Float, Column
from .connect import Base

__all__ = ['SextractorSchema']

class SextractorSchema(Base):
    """
    Hugs Pipeline SExtractor catalog schema.
    """

    __tablename__ = 'hugs_schema'
    
    id  = Column(Integer, primary_key=True)

    x_image = Column(Float, nullable=False)
    y_image = Column(Float, nullable=False)

    x_hsc = Column(Float, nullable=False)
    y_hsc = Column(Float, nullable=False)

"""
    ra = Column(Float, index=True)
    dec = Column(Float, nullable=False)

    mag_aper = Column(ARRAY(Float), nullable=False)
    mu_aper = Column(ARRAY(Float), nullable=False)

    a_image = Column(ARRAY(Float), nullable=False)
    b_image = Column(ARRAY(Float), nullable=False)
    theta_image = Column(ARRAY(Float), nullable=False)
    fwhm_image = Column(ARRAY(Float), nullable=False)

    kron_radius = Column(ARRAY(Float), nullable=False)
    flux_radius = Column(ARRAY(Float), nullable=False)
    mag_auto = Column(ARRAY(Float), nullable=False)

    flags = Column(ARRAY(Float), nullable=False)
"""

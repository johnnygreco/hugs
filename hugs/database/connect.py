from __future__ import division, print_function

import os
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session

Session = scoped_session(sessionmaker(autoflush=True, autocommit=False))
Base = declarative_base()

def connect(db_fn, overwrite=False):
    """
    Connect to database.
    """

    if overwrite and os.path.isfile(db_fn):
        assert 'safe' not in db_fn, 'cannot delete safe file'
        os.remove(db_fn)

    engine = create_engine('sqlite:///'+db_fn)
    Session.configure(bind=engine)
    Base.metadata.bind = engine
    Base.metadata.create_all(engine)

    return engine

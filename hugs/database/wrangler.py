"""
Data wrangler for the hugs database
"""

from .tables import Run, Tract, Patch, Source
from .connect import connect, Session

__all__ = ['start_session', 'Wrangler']


def start_session(db_fn):
    engine = connect(db_fn, False)
    session = Session()
    return session


class Wrangler(object):

    def __init__(self, session):
        pass

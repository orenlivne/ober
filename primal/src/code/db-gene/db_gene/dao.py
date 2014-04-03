'''
============================================================
Data Access Object (DAO) base class.

Created on November 2, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
from sqlalchemy.orm.session import sessionmaker
from sqlalchemy.engine import create_engine
    
class Dao(object):
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def __init__(self, url, **kwargs):
        '''Initialize from a SQLAlchemy database engine.'''
        try:
            self.engine = create_engine(url, **kwargs)
            self.Session = sessionmaker(bind=self.engine)
            self.online = True
        except ImportError  as e:
            self.online = False

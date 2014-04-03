'''
============================================================
Test the SQLAlchemy ORM framework basic functions. 

Created on October 30, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest
from numpy.ma.testutils import assert_equal
from sqlalchemy import create_engine
from tests.db.entities import Base, User
from sqlalchemy.orm.session import sessionmaker

class TestSqlAlchemy(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Create an in-memory database.'''
        self.engine = create_engine('sqlite:///:memory:', echo=True)
        Base.metadata.create_all(self.engine)
        
    def tearDown(self):
        '''Drop the database.'''
        pass
    
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_simple_select(self):
        '''Run a simple SELECT statement.'''
        assert_equal(self.engine.execute("select 1").scalar(), 1, 'Wrong select result')

    def test_user_entity(self):
        '''User class queries.'''
        # print repr(User.__table__) #@UndefinedVariable
        ed_user = User('ed', 'Ed Jones', 'edspassword')
        assert_equal(ed_user.id, None)
        assert_equal(ed_user.name, 'ed')
        assert_equal(ed_user.password, 'edspassword')
        
        Session = sessionmaker(bind=self.engine)
        session = Session()
        session.add(ed_user)
        
        our_user = session.query(User).filter_by(name='ed').first()
        assert_equal(our_user.id, 1)
        assert_equal(our_user.name, 'ed')
        assert_equal(our_user.password, 'edspassword')
        assert_equal(our_user, ed_user, 'User saved and loaded does not produce the orignial object')

        session.add_all([User('wendy', 'Wendy Williams', 'foobar'),
                         User('mary', 'Mary Contrary', 'xxg527'),
                         User('fred', 'Fred Flinstone', 'blah')])
        
        # Change the user object. Note that ed_user, our_user are the same reference - I guess you
        # can only bound one object to an entity identifier within a session
        ed_user.password = 'f8s7ccs'
        our_user.password = 'f8s7ccs1'
        assert_equal(len(session.dirty), 1, 'Ed user object changed but not in dirty object list')
        assert_equal(len(session.new), 3, 'New objects not in pending object list')
        
        # Committing changes saves the new password Ed set (the latest one set on the our_user reference)
        # and makes the ed_user object from transient to persistent 
        session.commit()
        assert_equal(ed_user.id, 1)
        our_user = session.query(User).filter_by(name='ed').first()
        assert_equal(our_user.password, 'f8s7ccs1')

        # Query        
        assert_equal([instance.name for instance in session.query(User).order_by(User.id)],
                     ['ed', 'wendy', 'mary', 'fred'])
        # Query for descriptors returns tuples (name,)
        assert_equal([name[0] for name in session.query(User.name).order_by(User.id)],
                     ['ed', 'wendy', 'mary', 'fred'])
        session.close()
        
    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

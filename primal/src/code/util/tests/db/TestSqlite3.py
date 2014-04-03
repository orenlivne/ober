'''
============================================================
Test the SQLITE3 database interaction library. 

Created on October 17, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import unittest, sqlite3
from numpy.ma.testutils import assert_equal

class TestSqlite3(unittest.TestCase):
    #---------------------------------------------
    # Constants
    #---------------------------------------------
    
    #---------------------------------------------
    # Setup Methods
    #---------------------------------------------
    def setUp(self):
        '''Create the database.'''
#        self.db_file = util.temp_filename(removeOnExit=False)
#        self.con = sqlite3.connect(self.db_file)        
        self.con = sqlite3.connect(":memory:")
        self.c = self.con.cursor()

        # Successful, con.commit() is called automatically afterwards
        with self.con:
            # Create table
            self.con.execute('''CREATE TABLE stocks
                             (date text, trans text, symbol text, qty real, price real)''')
            self.con.execute("create table person (id integer primary key, firstname varchar unique)")            

    def tearDown(self):
        '''Drop the database.'''
        # We can also close the conection if we are done with it.
        # Just be sure any changes have been committed or they will be lost.
        self.con.close()
#        os.remove(self.db_file)
        
    #---------------------------------------------
    # Test Methods
    #---------------------------------------------
    def test_insert_and_query(self):
        '''Check that a chain with printout filters works.'''
        (con, c) = (self.con, self.c)
    
        # Insert a row of data
        with con:
            con.execute("INSERT INTO stocks VALUES ('2006-01-05','BUY','RHAT',100,35.14)")
    
        # Do this instead
        t = ('RHAT',)
        c.execute('SELECT * FROM stocks WHERE symbol=?', t)
        assert_equal(c.fetchone(),
                     (u'2006-01-05', u'BUY', u'RHAT', 100.0, 35.14), 'Wrong single row query')
                
        # Larger example that inserts many records at a time
        purchases = [('2006-03-28', 'BUY', 'IBM', 1000, 45.00),
                     ('2006-04-05', 'BUY', 'MSFT', 1000, 72.00),
                     ('2006-04-06', 'SELL', 'IBM', 500, 53.00),
                    ]
        with con:
            con.executemany('INSERT INTO stocks VALUES (?,?,?,?,?)', purchases)
    
        assert_equal(list(con.execute('SELECT * FROM stocks ORDER BY price')),
                     [(u'2006-01-05', u'BUY', u'RHAT', 100.0, 35.14),
                      (u'2006-03-28', u'BUY', u'IBM', 1000.0, 45.0),
                      (u'2006-04-06', u'SELL', u'IBM', 500.0, 53.0),
                      (u'2006-04-05', u'BUY', u'MSFT', 1000.0, 72.0)], 'Wrong stocks query result')

    def test_con_as_context(self):
        '''Check that a chain with printout filters works.'''
        # Successful, con.commit() is called automatically afterwards
        con = self.con
        with con:
            con.execute("insert into person(firstname) values (?)", ("Joe",))
        
        # con.rollback() is called after the with block finishes with an exception, the
        # exception is still raised and must be caught
        try:
            with con:
                con.execute("insert into person(firstname) values (?)", ("Joe",))
            raise Exception('We shouldn''t here') 
        except sqlite3.IntegrityError:
            # Couldn't add Joe twice
            pass

    #---------------------------------------------
    # Private Methods
    #---------------------------------------------

'''
============================================================
Custom SQLAlchemy functions.

Created on October 30, 2012
@author: Oren Livne <livne@uchicago.edu>
@see http://docs.sqlalchemy.org/en/latest/orm/tutorial.html
============================================================
'''
from sqlalchemy.sql import expression
from sqlalchemy.ext.compiler import compiles
from sqlalchemy.types import Numeric

####################################################################################
'''GREATEST(a,b) - the maximum of a and b'''
class greatest(expression.FunctionElement):
    type = Numeric() #@ReservedAssignment
    name = 'greatest'

@compiles(greatest)
def default_greatest(element, compiler, **kw):
    return compiler.visit_function(element)

@compiles(greatest, 'sqlite')
@compiles(greatest, 'mssql')
@compiles(greatest, 'oracle')
def case_greatest(element, compiler, **kw):
    arg1, arg2 = list(element.clauses)
    return "CASE WHEN %s > %s THEN %s ELSE %s END" % (
        compiler.process(arg1),
        compiler.process(arg2),
        compiler.process(arg1),
        compiler.process(arg2),
    )

####################################################################################
'''LEAST(a,b) - the minimum of a and b'''
class least(expression.FunctionElement):
    type = Numeric() #@ReservedAssignment
    name = 'least'

@compiles(least)
def default_least(element, compiler, **kw):
    return compiler.visit_function(element)

@compiles(least, 'sqlite')
@compiles(least, 'mssql')
@compiles(least, 'oracle')
def case_least(element, compiler, **kw):
    arg1, arg2 = list(element.clauses)
    return "CASE WHEN %s < %s THEN %s ELSE %s END" % (
        compiler.process(arg1),
        compiler.process(arg2),
        compiler.process(arg1),
        compiler.process(arg2),
    )

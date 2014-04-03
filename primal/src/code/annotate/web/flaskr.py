#!/usr/bin/env python
'''
============================================================
Hutterites variant annotation web application.

Created on August 8, 2012
@author: Oren Livne <livne@uchicago.edu>
============================================================
'''
import os, sys, tempfile, subprocess as sub
from sqlite3 import dbapi2 as sqlite3
from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash
from werkzeug import secure_filename  # @UnresolvedImport
from flask.helpers import make_response

#---------------------------------------------
# Constants
#---------------------------------------------
ANNOTATION_EXEC='/home/oren/ober/system/bin-data/annotations/annotate-csv'

# create application
app = Flask(__name__)

# Load default config and override config from an environment variable
app.config.update(dict(
    DATABASE=os.path.join(app.root_path, 'flaskr.db'),
    DEBUG=True,
    SECRET_KEY='development key',
    USERNAME='admin',
    PASSWORD='default',
    UPLOAD_FOLDER=tempfile.gettempdir(),
    ALLOWED_EXTENSIONS=set(['csv'])
))
app.config.from_envvar('FLASKR_SETTINGS', silent=True)

#---------------------------------------------
# Database operations
#---------------------------------------------
def connect_db():
    '''Connects to the specific database.'''
    rv = sqlite3.connect(app.config['DATABASE'])
    rv.row_factory = sqlite3.Row
    return rv

def init_db():
    '''Creates the database tables.'''
    with app.app_context():
        db = get_db()
        with app.open_resource('schema.sql', mode='r') as f:
            db.cursor().executescript(f.read())
        db.commit()

def get_db():
    '''Opens a new database connection if there is none yet for the
    current application context.'''
    if not hasattr(g, 'sqlite_db'):
        g.sqlite_db = connect_db()
    return g.sqlite_db

@app.teardown_appcontext
def close_db(error):
    '''Closes the database again at the end of the request.'''
    if hasattr(g, 'sqlite_db'):
        g.sqlite_db.close()

def allowed_file(filename):
    '''is file''s extension allowed?'''
    return '.' in filename and filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

def annotated_file_name(filename):
    '''Transform an uploaded file name into the corresopnding downloaded file name.'''
    return os.path.basename(filename).rsplit('.', 1)[0] + '.annotated.csv'

def annotate_file(saved_file, chrom, bp):
    '''Business logic of annotated a file of Hutterite variants. Returns a return code,
    the output stream with of augmented CSV data and a list of errors.'''
    cmd = 'cat %s | %s %d %d' % (saved_file, ANNOTATION_EXEC, chrom, bp)
    return run_command(cmd)

def run_command(cmd, verbose=False):
    '''Run command in a sub-shell. Return the command's return code.'''
    try:
        p = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
        output, errors = p.communicate()
        if errors: sys.stderr.write('Errors: ' + repr(errors) + '\n')
        return p.returncode, output, filter(lambda x: x, errors.split('\n'))
    except ValueError:
        return -1, '', ['An unknown error annotating the file has occurred.']

#---------------------------------------------
# Routing Methods
#---------------------------------------------
#@app.route('/')
#def show_entries():
#    db = get_db()
#    cur = db.execute('select title, text from entries order by id desc')
#    entries = cur.fetchall()
#    return render_template('show_entries.html', entries=entries)

@app.route('/add', methods=['POST'])
def add_entry():
    if not session.get('logged_in'):
        abort(401)
    db = get_db()
    db.execute('insert into entries (title, text) values (?, ?)',
               [request.form['title'], request.form['text']])
    db.commit()
    flash('New entry was successfully posted')
    return redirect(url_for('show_entries'))

@app.route('/login', methods=['GET', 'POST'])
def login():
    error = None
    if request.method == 'POST':
        if request.form['username'] != app.config['USERNAME']:
            error = 'Invalid username'
        elif request.form['password'] != app.config['PASSWORD']:
            error = 'Invalid password'
        else:
            session['logged_in'] = True
            flash('Welcome, %s!' % (request.form['username'],))
            return redirect(url_for('annotate'))
    return render_template('login.html', error=error)

@app.route('/logout')
def logout():
    session.pop('logged_in', None)
    flash('You have been logged out.')
    return redirect(url_for('main_app'))

@app.route('/')
def main_app():
    return render_template('login.html')

@app.route('/annotate', methods=['GET', 'POST'])
def annotate():
    if not session.get('logged_in'): abort(401)
    errors = []
    if request.method == 'GET':
        return render_template('annotate.html', errors=errors)
    elif request.method == 'POST':
        # Form input validation
        try: chrom = int(request.form['chrom'])
        except ValueError: errors.append('Chromosome column number must be numeric.')
        
        try: bp = int(request.form['bp'])
        except ValueError: errors.append('Base-pair column number must be numeric.')
        
        uploaded_file = request.files['file']
        if not uploaded_file or not allowed_file(uploaded_file.filename):
            errors.append('Please upload a valid file ending with a .csv extension.')
        if errors: return render_template('annotate.html', errors=errors)

        # Upload file
        filename = secure_filename(uploaded_file.filename)
        saved_file = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        uploaded_file.save(saved_file)
#        flash('Saved file under ' + saved_file + '\n')
        returncode, output, run_errors = annotate_file(saved_file, chrom, bp)
        if returncode != 0 or run_errors:
            if not run_errors: errors.append('An unknown error annotating the file has occurred.')
            else:
                for error in run_errors: errors.append(error)            
            return render_template('annotate.html', errors=errors)
        else:
            response = make_response(output)
            response.headers['Content-Disposition'] = 'attachment; filename=%s' % (annotated_file_name(uploaded_file.filename),)
            os.remove(saved_file)
            return response

if __name__ == '__main__':
    init_db()
    app.run(debug=True)

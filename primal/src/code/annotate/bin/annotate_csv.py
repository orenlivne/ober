#!/usr/bin/env python
from flask import Flask, request

app = Flask(__name__)

@app.route('/')
def index():
    return "Ober Lab - Web Services"

@app.route('/annotate', methods=['POST', 'GET'])
def login():
    error = None
    if request.method == 'POST':
        if valid_login(request.form['username'],
                       request.form['password']):
            return log_the_user_in(request.form['username'])
        else:
            error = 'Invalid username/password'
    # the code below is executed if the request method
    # was GET or the credentials were invalid
    return render_template('login.html', error=error)
if __name__ == '__main__':
    app.run(host='0.0.0.0', debug = True)

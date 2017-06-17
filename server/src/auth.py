#!/usr/bin/env python
# -*- Mode: Python; tab-width: 4; indent-tabs-mode: nil; coding: utf-8; -*-
# vim:set ft=python ts=4 sw=4 sts=4 autoindent:

'''
Authentication and authorization mechanisms.

Author:     Pontus Stenetorp    <pontus is s u-tokyo ac jp>
            Illes Solt          <solt tmit bme hu>
Version:    2011-04-21
'''

from hashlib import sha512
from os.path import dirname, join as path_join, isdir
from os import mkdir
import mysql.connector
from sys import stderr
from config import DATA_DIR

try:
    from os.path import relpath
except ImportError:
    # relpath new to python 2.6; use our implementation if not found
    from common import relpath
from common import ProtocolError
from config import USER_PASSWORD, DATA_DIR
from message import Messager
from session import get_session, invalidate_session
from projectconfig import ProjectConfiguration

mysql_config = {
  'user': 'oshin',
  'password': 'bratmedline',
  'host': 'east-mysql-instance1.ckwuvvjrdeqm.us-east-1.rds.amazonaws.com',
  'port': 3306,
  'database': 'medline',
}

# To raise if the authority to carry out an operation is lacking
class NotAuthorisedError(ProtocolError):
    def __init__(self, attempted_action):
        self.attempted_action = attempted_action

    def __str__(self):
        return 'Login required to perform "%s"' % self.attempted_action

    def json(self, json_dic):
        json_dic['exception'] = 'notAuthorised'
        return json_dic


# File/data access denial
class AccessDeniedError(ProtocolError):
    def __init__(self):
        pass

    def __str__(self):
        return 'Access Denied'

    def json(self, json_dic):
        json_dic['exception'] = 'accessDenied'
        # TODO: Client should be responsible here
        Messager.error('Access Denied')
        return json_dic


class InvalidAuthError(ProtocolError):
    def __init__(self):
        pass

    def __str__(self):
        return 'Incorrect login and/or password'

    def json(self, json_dic):
        json_dic['exception'] = 'invalidAuth'
        return json_dic


def _is_authenticated(user, password):

    if not user or not password:
        return False

    cnx = mysql.connector.connect(**mysql_config);
    cursor = cnx.cursor()
    query = "select password from users where username='%s'" % user;
    cursor.execute(query)
    correct_password = ""
    for val in cursor:
       correct_password = val[0]
    cursor.close()
    cnx.close()
    return password==correct_password

def _password_hash(password):
    return sha512(password).hexdigest()

def login(user, password):
    if not _is_authenticated(user, password):
        raise InvalidAuthError

    get_session()['user'] = user
    Messager.info('Hello!')
    return {}

def logout():
    try:
        del get_session()['user']
    except KeyError:
        # Already deleted, let it slide
        pass
    # TODO: Really send this message?
    Messager.info('Bye!')
    return {}

def whoami():
    json_dic = {}
    try:
        json_dic['user'] = get_session().get('user')
    except KeyError:
        # TODO: Really send this message?
        Messager.error('Not logged in!', duration=3)
    return json_dic

def allowed_to_read(real_path):
    data_path = path_join('/', relpath(real_path, DATA_DIR))
    # add trailing slash to directories, required to comply to robots.txt
    if isdir(real_path):
        data_path = '%s/' % ( data_path )
        
    real_dir = dirname(real_path)
    robotparser = ProjectConfiguration(real_dir).get_access_control()
    if robotparser is None:
        return True # default allow

    try:
        user = get_session().get('user')
    except KeyError:
        user = None

    if user is None:
        user = 'guest'

    #display_message('Path: %s, dir: %s, user: %s, ' % (data_path, real_dir, user), type='error', duration=-1)

    return robotparser.can_fetch(user, data_path)

def signup(user, password):

    if not user or not password:
        Messager.info('Username or password not provided!')
        return {}

    cnx = mysql.connector.connect(**mysql_config);
    cursor = cnx.cursor()
    query = "select username from users where username='%s'" % user;
    cursor.execute(query)
    numrows = 0 
    for val in cursor:
       numrows += 1

    if numrows>0: 
       cursor.close()
       cnx.close()
       Messager.info('User already exists!')
       return {}

    query = "insert into users(username,password) values('%s','%s')" % (user,password);
    cursor.execute(query)
    cursor.close()
    cnx.commit()
    cnx.close()
    mkdir( DATA_DIR+'/'+user, 0755 );
    mkdir( DATA_DIR+'/'+user+'/old', 0755 );
    f = open(DATA_DIR+'/'+user+'/annotation.conf','w')
    f.write("[entities]\n\nPatient\nIntervention\nOutcome\n\n[relations]\n[events]\n[attributes]\n")
    f.close()
    Messager.info('Signed up! Hello!')
    return {}

# TODO: Unittesting

###############################################################################
## Stigmee: A 3D browser and decentralized social network.
## Copyright 2021 Duron Alain <duron.alain@gmail.com>
##
## This file is part of Stigmee.
##
## Project : Stigmee BeeBot
## Version : 0.0-1
## Date : 20-11-2021
## Author : Alain Duron
## File : BeeDB.py
##
## Stigmee is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

# TODO:
#   Handle other database systems if necessary

# Imports
try:
    import mariadb
    import sys
except ImportError:
    print("[E] can't find the 'google' module. try : 'pip install mariadb'")


"""
BeeDB is a class that encapsulates MariaDB connector, simplifying Database connection
for other Stigmee BeeBot modules
"""
class BeeDB:

    def __init__(self, _user, _password, _host, _port, _database):

        # Connect to MariaDB Platform
        self.user = _user
        self.password = _password
        self.host = _host
        self.port = _port
        self.database = _database
        self.cur = None
        self.conn = None

        # Setup the connection
        self.connect()
        self.get_cursor()

    def connect(self):

        try:
            self.conn = mariadb.connect(
                user=self.user,
                password=self.password,
                host=self.host,
                port=self.port,
                database=self.database
            )
        except mariadb.Error as e:
            print(f"Error connecting to MariaDB Platform: {e}")
            sys.exit(1)

        # By default auto commit is on
        # disable it and make it an argument in execute_stmt
        self.conn.autocommit = False

    def get_cursor(self):

        # Get Cursor
        self.cur = self.conn.cursor()

    def execute_query(self, _query, _paramtuple=None):
        # example :
        # _query : "SELECT first_name,last_name FROM employees WHERE first_name=?",
        # _paramtuple : (some_name,)

        if _paramtuple is not None:
            self.cur.execute(_query,_paramtuple)
        else:
            self.cur.execute(_query)

        return self.cur

    def execute_stmt(self, _stmt, _paramtuple=None, _commit=False):

        if _paramtuple is not None:
            self.cur.execute(_stmt,_paramtuple)
        else:
            self.cur.execute(_stmt)

        if _commit:
            self.commit()

    def commit(self):

        self.conn.commit()

    def rollback(self):

        self.conn.rollback()

    def close(self):

        self.conn.close()
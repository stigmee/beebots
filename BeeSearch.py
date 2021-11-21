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
## File : BeeSearch.py
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
#   Add other search engine APIs for input
#   Provide better integration with other classes

# Imports
try:
    from googlesearch import search
except ImportError:
    print("[E] can't find the 'google' module. try : 'pip install google'")


class GoogleSearch:

    def __init__(self, _query, _tld='com', _lang='fr', _num=200, _start=0, _stop=200, _verbose=False):
        """
        :param _query: The keyword(s) to be search for, eg "Economie Bleue"
        :param _tld: the top level domain where the query should be executed (.com, .fr, .co.in...)
        :param _lang: the query language (en, fr...)
        :param _num: the number of results to retreive
        :param _start: the starting index for retreived results
        :param _stop: the stoping index
        :param _verbose: prints

        The pause option of the API will force the program to wait between HTTP requests.
        If set too low, google may block the caller's IP, so better keep that to at least 2.
        """

        # Inputs
        self.query = _query
        self.verbose = _verbose
        self.tld = _tld
        self.lang = _lang
        self.num = _num
        self.start = _start
        self.stop = _stop

        # Output structure

        self.results = []

        # Execute the query and initialize an indexed result structure

        self.runSearch()

    def getResults(self, _query=None, _force=False):

        if (_query is None) and _force:
            # re run the search anyway
            self.runSearch()
            return self.results
        elif (_query is None) or (self.query == _query and _force is False):
            # By default I return the existing cache, if the query is called using the same keywords
            # this avoids doing multiple times the same query, which can also trigger blocking
            return self.results
        else:
            # Update the query and search
            self.query = _query
            self.runSearch()
            return self.results

    def runSearch(self):
        for i, url in enumerate(search(self.query, tld=self.tld, lang=self.lang, num=self.num,
                                       start=self.start, stop=self.stop, pause=2)):
            if self.verbose: print("[r] {} : {}".format(i + self.start, url))
            self.results.append((i + self.start, url))

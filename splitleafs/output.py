#!/usr/bin/env python

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""
Implement the split algorithms.
"""

from __future__ import print_function, division, with_statement
import textwrap
import sys

__author__ = "Jonathan Barnoud"


def write_ndx(groups):
    """
    Write a gromacs index file with the given groups.
    """
    for group_name, atomids in groups.items():
        print("[ {0} ]".format(group_name))
        group_str = " ".join([str(i) for i in atomids])
        print("\n".join(textwrap.wrap(group_str, 80)))


def stats(groups, outfile=sys.stderr):
    """
    Display some statistics on the leaflets.
    """
    for group_name, atomids in groups.items():
        print("{0}: {1} atoms".format(group_name, len(atomids)),
              file=outfile)
    if len(groups["upper_leaflet"]) == len(groups["lower_leaflet"]):
        print("The membrane is symmetric.", file=outfile)

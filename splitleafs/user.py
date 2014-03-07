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
Handle user inputs.
"""

from __future__ import print_function, division, with_statement
import argparse
import os
from .selection import parse_selection

__author__ = "Jonathan Barnoud"


def isfile(path):
    """
    Check if path is an existing file.

    If not, raise an error. Else, return the path.
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_options(argv):
    """
    Read the command line arguments.
    """
    usage = ("%(prog)s [options] < input > output.ndx\n"
             "       %(prog)s [options] -- input > output.ndx")
    parser = argparse.ArgumentParser(description=__doc__, usage=usage)
    parser.add_argument("input", default=None, nargs='?', type=isfile,
                        help="The input structure.")
    parser.add_argument("--axis", "-d", choices="xyz", default="z",
                        help="Axis normal to the bilayer.")
    parser.add_argument("--atom", "-a", type=parse_selection,
                        default=[("P1",)], nargs='+',
                        help="Reference atom name.")
    parser.add_argument("--format", "-f", type=str,
                        default="auto", choices=["gro", "pdb", "auto"],
                        help="Input file format.")
    keep_options = parser.add_mutually_exclusive_group()
    keep_options.add_argument("--keep-residue", "-r", action="store_true",
                              dest="keep_residue", default=False,
                              help="Keep the whole residues.")
    keep_options.add_argument("--keep-atom", "-k", action="store_false",
                              dest="keep_residue", default=False,
                              help="Keep only the atom of reference.")
    args = parser.parse_args(argv)
    return args



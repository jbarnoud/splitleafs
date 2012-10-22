#!usr/bin/env python

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
Test suite for splitleafs.
"""

from unittest import TestCase, main
import os
import splitleafs

__author__ = "Jonathan Barnoud"

REFDIR = "test_resources"
PRECISION = 4


class TestLibrary(TestCase):
    """
    Test the functions.
    """
    def test_read_gro(self):
        """
        Test gro file reading on a very simple case.
        """
        path = os.path.join(REFDIR, "simple.gro")
        reference = [
            {"resid":1, "resname":"FOO", "atom_name":"BAR", "atomid":1,
             "x":0.0, "y":1.0, "z":2.0},
            {"resid":1, "resname":"FOO", "atom_name":"BAR2", "atomid":2,
             "x":3.0, "y":4.0, "z":5.0},
            {"resid":2, "resname":"OOF", "atom_name":"BAR3", "atomid":3,
             "x":6.0, "y":7.0, "z":-8.0},
        ]
        with open(path) as infile:
            atoms = splitleafs.read_gro(infile.readlines()[2:-1])
        for ref, atom in zip(reference, atoms):
            for key, value in ref.iteritems():
                if not key in "xyz":
                    self.assertEqual(value, atom[key],
                                     ("Different value for key {0}: {1} "
                                      "instead of {2}")
                                     .format(key, atom[key], value))
                else:
                    self.assertAlmostEqual(value, atom[key], PRECISION,
                                           ("Different value for key {0}: {1} "
                                            "instead of {2}")
                                           .format(key, atom[key], value))


if __name__ == "__main__":
    main()

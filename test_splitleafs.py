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

from __future__ import print_function
from unittest import TestCase, main
import os
import splitleafs
import subprocess

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
            for key, value in ref.items():
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

    def test_split_get_res(self):
        """
        Test if the group formed by the split_get_res function are consistent.
        """
        # Split the bilayer
        path = os.path.join(REFDIR, "membrane.gro")
        axis = "z"
        atoms = list(splitleafs.read_gro(open(path).readlines()[2:-1]))
        selection = list(splitleafs.select_atom_name(atoms, "P1"))
        coordinates = splitleafs.axis_coordinates(selection, axis)
        average = splitleafs.mean(coordinates)
        groups = splitleafs.split_get_res(atoms, average, axis, "P1")

        # Do the group have the right size?
        natoms_popc = sum((1 for atom in atoms if atom["resname"] == "POPC"))
        print("natom_popc: {0}".format(natoms_popc))
        for key, value in groups.items():
            self.assertEqual(len(value) * 2, natoms_popc,
                             ("The {0} group contains {1} atoms it should "
                              "contain {2} ones.")
                             .format(key, len(value), natoms_popc // 2))

        # Is there doubles in the groups?
        for key, value in groups.items():
            self.assertEqual(len(value), len(set(value)),
                             "There is double values in the {0} group."
                             .format(key))

        # Is the group sorted? A group can not be sorted if the input file
        # is not but the test input file is correctly sorted.
        for key, value in groups.items():
            ordered = value[:]
            ordered.sort()
            self.assertEqual(value, ordered,
                             "The {0} group is nor correctly ordered."
                             .format(key))

        # Are the groups the same as the reference?
        ndx_reference = os.path.join(REFDIR, "leafs.ndx")
        with open(ndx_reference) as infile:
            reference = read_ndx(infile)
        for key, value in reference.items():
            self.assertEqual(value, groups[key],
                             ("The {0} group is different between the "
                              "function output and the reference.")
                             .format(key))


class TestProgram(TestCase):
    """
    Test that the program run correctly.
    """
    def test_run(self):
        """
        Launch the program.
        """
        with open("test_output.ndx", "w") as out:
            with open("test_error.txt", "w") as err:
                status = subprocess.call(["./splitleafs.py",
                                          "test_resources/membrane.gro", "-r"],
                                         stdout=out, stderr=err)
        with open("test_error.txt") as err:
            for line in err:
                print(err)
        self.assertEqual(status, 0, "Error while calling the program.")
        ndx_reference = os.path.join(REFDIR, "leafs.ndx")
        with open(ndx_reference) as infile:
            reference = read_ndx(infile)
        with open("test_output.ndx") as infile:
            groups = read_ndx(infile)
        self.assertEqual(reference, groups, "Program output is wrong.")

    def test_fail_run(self):
        """
        Test that the test suit actually catch crashes.
        """
        with open("test_output.ndx", "w") as out:
            with open("test_error.txt", "w") as err:
                status = subprocess.call(["./splitleafs.py",
                                          "test_resources/non-existing.gro",
                                          "-r"], stdout=out, stderr=err)
        with open("test_error.txt") as err:
            for line in err:
                print(err)
        self.assertNotEqual(status, 0, "The program should have crash.")


def read_ndx(infile):
    """
    Read a GROMACS index file and return a dictionary.

    :Parameters:
        - infile : a file descriptor like instance of a ndx file

    :Return:
        - A dictionary like {group name : list of indices}.
    """
    indices = {}
    current_group = None
    groups = []
    for line in infile:
        # Remove comments if any
        comment_start = line.find(";")
        if comment_start > -1:
            line = line[:comment_start]
        if "[" in line:
            current_group = line
            current_group = current_group.replace("[", "")
            current_group = current_group.replace("]", "")
            current_group = current_group.strip()
            groups.append(current_group)
            indices[current_group] = []
        elif not current_group is None:
            indices[current_group] += [int(i) for i in line.split()]
    return indices


if __name__ == "__main__":
    main()

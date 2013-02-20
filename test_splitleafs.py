#!usr/bin/env python

#    This file is part of splitleafs.
#
#    Splitleafs is free software: you can redistribute it and/or modify
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
import contextlib
import itertools
import os
import splitleafs
import subprocess
import sys

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
            atoms = list(splitleafs.read_gro(infile))
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
        with open(path) as infile:
            atoms = list(splitleafs.read_gro(infile))
        selection = list(splitleafs.select_atoms(atoms, [("P1",)]))
        coordinates = splitleafs.axis_coordinates(selection, axis)
        average = splitleafs.mean(coordinates)
        groups = splitleafs.split_get_res(atoms, average, axis, [("P1",)])

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

    def test_multiple_references(self):
        """
        Test if it works with several reference atoms in the same residue.
        """
        # Split the bilayer using POPC:P1 and POPC:C12 as reference atoms
        path = os.path.join(REFDIR, "membrane.gro")
        axis = "z"
        select_code = [("POPC", "P1"), ("POPC", "C12")]
        with open(path) as infile:
            atoms = list(splitleafs.read_gro(infile))
        print(len(atoms))
        selection = list(splitleafs.select_atoms(atoms, select_code))
        print(len(selection))
        coordinates = splitleafs.axis_coordinates(selection, axis)
        average = splitleafs.mean(coordinates)
        groups = splitleafs.split_get_res(atoms, average, axis, select_code)

        # The resulting groups should be the exact same as the one from the
        # reference output built with just POPC:P1 as reference atom
        ndx_reference = os.path.join(REFDIR, "leafs.ndx")
        with open(ndx_reference) as infile:
            reference = read_ndx(infile)
        for key, value in reference.items():
            self.assertEqual(value, groups[key],
                             ("The {0} group is different between the "
                              "function output and the reference.")
                             .format(key))



    def test_keep_options(self):
        """
        Test that the program complains if both -r and -k are used.
        """
        with _redirect_stderr(sys.stdout):
            # With both -r and -k the function should crash
            argv = ["-r", "-k", "{0}/membrane.gro".format(REFDIR)]
            self.assertRaises(SystemExit, splitleafs.get_options, *[argv])
            # With only -r or -k it should not crash
            argv = ["-r", "{0}/membrane.gro".format(REFDIR)]
            splitleafs.get_options(argv)
            argv = ["-k", "{0}/membrane.gro".format(REFDIR)]
            splitleafs.get_options(argv)


class TestProgram(TestCase):
    """
    Test that the program run correctly.
    """
    def test_run_residue(self):
        """
        Launch the program and keep residues.
        """
        # This nested with statements are kept for compatibility with
        # python 2.6 otherwise it could be written
        # with open("test_output.ndx", "w") as out, \
        #      open("test_error.txt", "w") as err:
        with open("test_output.ndx", "w") as out:
            with open("test_error.txt", "w") as err:
                status = subprocess.call(["./splitleafs.py",
                                          "test_resources/membrane.gro", "-r"],
                                         stdout=out, stderr=err)
        with open("test_error.txt") as err:
            for line in err:
                print(line)
        self.assertEqual(status, 0, "Error while calling the program.")
        ndx_reference = os.path.join(REFDIR, "leafs.ndx")
        with open(ndx_reference) as infile:
            reference = read_ndx(infile)
        with open("test_output.ndx") as infile:
            groups = read_ndx(infile)
        self.assertEqual(reference, groups, "Program output is wrong.")

    def test_run_no_residue(self):
        """
        Launch the program and do not keep residues.
        """
        # This nested with statements are kept for compatibility with
        # python 2.6
        with open("test_output.ndx", "w") as out:
            with open("test_error.txt", "w") as err:
                status = subprocess.call(["./splitleafs.py",
                                          "test_resources/membrane.gro"],
                                         stdout=out, stderr=err)
        with open("test_error.txt") as err:
            for line in err:
                print(line)
        self.assertEqual(status, 0, "Error while calling the program.")
        ndx_reference = os.path.join(REFDIR, "leafsP1.ndx")
        with open(ndx_reference) as infile:
            reference = read_ndx(infile)
        with open("test_output.ndx") as infile:
            groups = read_ndx(infile)
        self.assertEqual(reference, groups, "Program output is wrong.")

    def test_fail_run(self):
        """
        Test that the test suit actually catch crashes.
        """
        # This nested with statements are kept for compatibility with
        # python 2.6
        with open("test_output.ndx", "w") as out:
            with open("test_error.txt", "w") as err:
                status = subprocess.call(["./splitleafs.py",
                                          "test_resources/non-existing.gro",
                                          "-r"], stdout=out, stderr=err)
        with open("test_error.txt") as err:
            for line in err:
                print(line)
        self.assertNotEqual(status, 0, "The program should have crash.")

    def test_begin_extra_lines(self):
        """
        There should not be anything befole the first group header.
        """
        # This nested with statements are kept for compatibility with
        # python 2.6
        with open("test_output.ndx", "w") as out:
            with open("test_error.txt", "w") as err:
                subprocess.call(["./splitleafs.py",
                                 "test_resources/membrane.gro"],
                                stdout=out, stderr=err)
        with open("test_output.ndx") as exec_out:
            extra_lines = 0
            for _ in itertools.takewhile(lambda x: "[" not in x, exec_out):
                extra_lines += 1
        self.assertEqual(extra_lines, 0,
                         ("There is {0} extra lines at the beginning of the "
                          "standard output.".format(extra_lines)))


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


@contextlib.contextmanager
def _redirect_stderr(destination=sys.stdout):
    """
    Redirect sys.stderr to an other file descriptor (sys.stdout by default).

    This function is a context manager and is used like as followed::

        >>> with _redirect_stderr():
        ...     print >> sys.stderr, "Something"

    In this example the print is done in sys.stdout even if sys.stderr is
    explicitly mentioned.

    Stderr can also be redirected to a file descriptor::

        >>> with open("my_file", "wt") as outfile:
        ...     with _redirect_stderr(outfile):
        ...         print >> sys.stderr, "Something"

    """
    old_stderr = sys.stderr
    sys.stderr = destination
    yield
    sys.stderr = old_stderr


if __name__ == "__main__":
    main()

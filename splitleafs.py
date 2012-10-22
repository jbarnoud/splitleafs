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
Write a GROMACS index file with one group per membrane leaflet.
"""

from __future__ import print_function, division
import argparse
import textwrap
import sys

__author__ = "Jonathan Barnoud"

GRO_FIELDS = {
    "resid": ((0, 5), int),
    "resname": ((5, 10), str),
    "atom_name": ((10, 15), str),
    "atomid": ((15, 20), int),
    "x": ((20, 28), float),
    "y": ((28, 36), float),
    "z": ((36, 44), float),
}


def read_gro(lines):
    """
    Read the atoms from a gro file.

    This function create a generator that yield a dictionary per line of the
    gro file.
    """
    for line in lines:
        atom = dict(((key, convert(line[begin:end].strip()))
                    for key, ((begin, end), convert)
                    in GRO_FIELDS.iteritems()))
        yield atom


def select_atom_name(atoms, atom_name):
    """
    Select only the atoms with the given atom name.

    This function create a generator that yield the dictionary for the atoms
    with the given atom name.
    """
    for atom in atoms:
        if atom["atom_name"] == atom_name:
            yield atom


def axis_coordinates(atoms, axis):
    """
    Get the coordinate of the atom along the given axis.

    Create a generator on the atom coordinate along the axis of interest.
    """
    for atom in atoms:
        yield atom[axis]


def mean(values):
    """
    Calculate the mean of an iterator.
    """
    summation = 0
    nelements = 0
    for value in values:
        summation += value
        nelements += 1
    return summation / nelements


def split(atoms, average, axis):
    """
    Split the leaflets along the given axis.
    """
    groups = {"upper_leaflet": [], "lower_leaflet": []}
    for atom in atoms:
        if atom[axis] >= average:
            groups["upper_leaflet"].append(atom["atomid"])
        else:
            groups["lower_leaflet"].append(atom["atomid"])
    return groups


def split_get_res(atoms, average, axis, atom_name):
    groups = {"upper_leaflet": [], "lower_leaflet": []}
    current_res = None
    current_group = None
    for atom in atoms:
        if atom["atom_name"] == atom_name:
            current_res = atom["resid"]
            if atom[axis] >= average:
                current_group = "upper_leaflet"
            else:
                current_group = "lower_leaflet"
            groups[current_group].append(atom["atomid"])
        if not current_res is None and atom["resid"] == current_res:
            groups[current_group].append(atom["atomid"])
    return groups


def write_ndx(groups):
    """
    Write a gromacs index file with the given groups.
    """
    for group_name, atomids in groups.iteritems():
        print("[ {0} ]".format(group_name))
        group_str = " ".join([str(i) for i in atomids])
        print("\n".join(textwrap.wrap(group_str, 80)))


def split_leaflets(lines, axis, atom_name, res=False):
    """
    Split bilayer leaflets from a gromacs gro file along the given axis.
    """
    axis = axis.lower()
    atoms = list(read_gro(lines))
    selection = list(select_atom_name(atoms, atom_name))
    coordinates = axis_coordinates(selection, axis)
    average = mean(coordinates)
    if res:
        groups = split_get_res(atoms, average, axis, atom_name)
    else:
        groups = split(selection, average, axis)
    write_ndx(groups)
    return groups


def get_options(argv):
    """
    Read the command line arguments.
    """
    usage = "%(prog)s [options] < input.gro > output.ndx"
    parser = argparse.ArgumentParser(description=__doc__, usage=usage)
    parser.add_argument("--axis", "-d", choices="xyz", default="z",
                        help="Axis normal to the bilayer.")
    parser.add_argument("--atom", "-a", type=str, default="P1",
                        help="Reference atom name.")
    parser.add_argument("--keep_residue", "-r", action="store_true",
                        dest="keep_residue", default=False,
                        help="Keep the whole residues.")
    parser.add_argument("--keep_atom", "-k", action="store_false",
                        dest="keep_residue", default=False,
                        help="Keep only the atom of reference.")
    args = parser.parse_args(argv)
    return args


def main():
    """
    Run every thing from the command line.
    """
    args = get_options(sys.argv[1:])
    lines = sys.stdin.readlines()[2:-1]
    groups = split_leaflets(lines, args.axis, args.atom, args.keep_residue)
    for group_name, atomids in groups.iteritems():
        print("{0}: {1} atoms".format(group_name, len(atomids)),
              file=sys.stderr)

if __name__ == "__main__":
    main()

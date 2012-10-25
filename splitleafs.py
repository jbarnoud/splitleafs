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
import itertools
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

PDB_FIELDS = {
    "resid": ((22, 26), int),
    "resname": ((18, 20), str),
    "atom_name": ((12, 16), str),
    "atomid": ((6, 11), int),
    "x": ((30, 38), float),
    "y": ((38, 46), float),
    "z": ((46, 54), float),
}

PDB_SECTIONS = (
    "HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND", "SOURCE",
    "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE",
    "JRNL", "REMARK", "DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES",
    "MODRES", "HET", "FORMUL", "HETNAM", "HETSYN", "HELIX", "SHEET", "SSBOND",
    "LINK", "CISPEP", "SITE", "CRYST1", "MTRIX", "ORIGX", "SCALE", "MODEL",
    "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL", "CONECT", "MASTER", "END",
)


class FormatError(Exception):
    """
    Exception raised when the file format is wrong.
    """
    pass


def read_gro(lines):
    """
    Read the atoms from a gro file.

    This function create a generator that yield a dictionary per line of the
    gro file.

    :Parameters:
        - lines: an iterator over atom lines from the gro file. The two header
                 lines and the bottom line describing the box have to be
                 excluded.

    :Raise:
        - FormatError: raised if the file format does not fit.
    """
    for line in lines:
        try:
            atom = dict(((key, convert(line[begin:end].strip()))
                        for key, ((begin, end), convert)
                        in GRO_FIELDS.iteritems()))
        except ValueError:
            raise FormatError
        yield atom


def read_pdb(lines):
    """
    Read the atoms from a pdb file.

    This function create a generator that yield a dictionary per line of the
    pdb file.

    :Parameters:
        - lines: an iterator over the lines from the PDB file

    :Raise:
        - FormatError: raised if the file format does not fit.
    """
    for line in lines:
        # Test if the line starts as it should in a PDB file
        if not (line[0:6].strip() in PDB_SECTIONS or line.strip() == ""):
            # MTRIX, ORIGX and SCALE pdb sections can be followed be a random
            # number so they will trigger the previous test
            if not line[0:4] in ("MTRIX", "ORIGX", "SCALE"):
                raise FormatError('PDB line should not start with "{0}"'
                                  .format(line[0:6]))
        if line[0:6] == "ATOM  ":
            try:
                atom = dict(((key, convert(line[begin:end].strip()))
                            for key, ((begin, end), convert)
                            in PDB_FIELDS.iteritems()))
            except ValueError:
                raise FormatError
            yield atom


def select_atom_name(atoms, atom_name):
    """
    Select only the atoms with the given atom name.

    This function create a generator that yield the dictionary for the atoms
    with the given atom name.

    :Parameters:
        - atoms: an iterator over the atom dictionaries
        - atom_name: the name of the atoms to keep
    """
    for atom in atoms:
        if atom["atom_name"] == atom_name:
            yield atom


def axis_coordinates(atoms, axis):
    """
    Get the coordinate of the atom along the given axis.

    Create a generator on the atom coordinate along the axis of interest.

    :Parameters:
        - atoms:  an iterator over the atom dictionaries
        - axis: the name of the dimension normal to the membrane (x, y or z)
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
    """
    Split the leaflets along the given axis and keep the whole residue.
    """
    groups = {"upper_leaflet": [], "lower_leaflet": []}
    keep_res = None
    current_resid = None
    current_res_atoms = []
    current_group = None
    for atom in atoms:
        # Keep track of the atoms of the residue, this is needed to have
        # have in the groups the atoms of a residue of interest that have been
        # read before the reference atom
        if atom["resid"] != current_resid:
            # We start a new residue
            current_resid = atom["resid"]
            current_res_atoms = []
            # Always reset the keep_res variable when we change the residue.
            # If we miss that we can catch extra residues because of the
            # reseting of the residue number that happend when the resid become
            # too big.
            keep_res = None
        current_res_atoms.append(atom["atomid"])
        # Split the residues of interest
        if atom["atom_name"] == atom_name:
            keep_res = atom["resid"]
            # Choose the group
            if atom[axis] >= average:
                current_group = "upper_leaflet"
            else:
                current_group = "lower_leaflet"
            # Add the atom of the residue that were read before the reference
            # atom. Do not include the last atom of the list since it is the
            # current atom and that he will be added to the group later.
            groups[current_group] += current_res_atoms[:-1]
        # Store the atom in the right group
        if not keep_res is None and atom["resid"] == keep_res:
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


def split_leaflets(infile, axis, atom_name, res=False, input_format="gro"):
    """
    Split bilayer leaflets from a gromacs gro file along the given axis.
    """
    axis = axis.lower()
    if input_format == "gro":
        atoms = list(read_gro(list(infile)[2:-1]))
    else:
        atoms = list(read_pdb(infile))
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
    usage = ("%(prog)s [options] < input > output.ndx\n"
             "       %(prog)s [options] input > output.ndx")
    parser = argparse.ArgumentParser(description=__doc__, usage=usage)
    parser.add_argument("input", default=None, nargs='?',
                        help="The input structure.")
    parser.add_argument("--axis", "-d", choices="xyz", default="z",
                        help="Axis normal to the bilayer.")
    parser.add_argument("--atom", "-a", type=str, default="P1",
                        help="Reference atom name.")
    parser.add_argument("--format", "-f", type=str,
                        default="auto", choices=["gro", "pdb", "auto"],
                        help="Input file format.")
    parser.add_argument("--keep_residue", "-r", action="store_true",
                        dest="keep_residue", default=False,
                        help="Keep the whole residues.")
    parser.add_argument("--keep_atom", "-k", action="store_false",
                        dest="keep_residue", default=False,
                        help="Keep only the atom of reference.")
    args = parser.parse_args(argv)
    return args


def guess_format(infile):
    """
    Guess the format of the input file among gro and pdb.

    Look if the file is a PDB one or assume it is a gro file.

    Return the format and an iterator that mimic the input file.
    """
    first_line = infile.readline()
    if (first_line[0:6].strip() in PDB_SECTIONS
            or first_line[0:5] in ("MTRIX", "ORIGX", "SCALE")):
        input_format = "pdb"
    else:
        input_format = "gro"
    mod_infile = itertools.chain([first_line], infile)
    return input_format, mod_infile


def main():
    """
    Run everything from the command line.
    """
    args = get_options(sys.argv[1:])
    if args.input is None:
        infile = sys.stdin
    else:
        infile = open(args.input)
    with infile:
        # Guess the format
        if args.format == "auto":
            input_format, mod_infile = guess_format(infile)
        else:
            input_format = args.format
            mod_infile = infile
        # Do the work
        try:
            groups = split_leaflets(mod_infile, args.axis, args.atom,
                                    args.keep_residue, input_format)
        # Complain if the format is wrong
        except FormatError:
            print(("Error while reading the input. Are you sure your file is "
                   "in the {0} format?").format(args.format), file=sys.stderr)
            sys.exit(1)
        # Complain if the reference atom is absent
        except ZeroDivisionError:
            print(("The reference atom looks abent from your input. Are you "
                   "sure there is some {0} atoms in your system?")
                  .format(args.atom), file=sys.stderr)
            sys.exit(1)
    # Display the number of atoms per group
    for group_name, atomids in groups.iteritems():
        print("{0}: {1} atoms".format(group_name, len(atomids)),
              file=sys.stderr)

if __name__ == "__main__":
    main()

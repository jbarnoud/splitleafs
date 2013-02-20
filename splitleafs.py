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
Write a GROMACS index file with one group per membrane leaflet.
"""

from __future__ import print_function, division, with_statement
import argparse
import itertools
import textwrap
import sys
import os

__author__ = "Jonathan Barnoud"

# File format description. The key is the name of the field, the value is a
# tuple from which the first element is the first (included) and last
# (excluded) indices of the field in the line, and the second element the type
# of the field content.
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
    "resname": ((17, 20), str),
    "atom_name": ((12, 16), str),
    "atomid": ((6, 11), int),
    "x": ((30, 38), float),
    "y": ((38, 46), float),
    "z": ((46, 54), float),
}

# All the authorized values for the first field of a line in a PDB file
PDB_SECTIONS = (
    "HEADER", "OBSLTE", "TITLE", "SPLT", "CAVEAT", "COMPND", "SOURCE",
    "KEYWDS", "EXPDTA", "NUMMDL", "MDLTYP", "AUTHOR", "REVDAT", "SPRSDE",
    "JRNL", "REMARK", "DBREF", "DBREF1", "DBREF2", "SEQADV", "SEQRES",
    "MODRES", "HET", "FORMUL", "HETNAM", "HETSYN", "HELIX", "SHEET", "SSBOND",
    "LINK", "CISPEP", "SITE", "CRYST1", "MTRIX", "ORIGX", "SCALE", "MODEL",
    "ATOM", "ANISOU", "TER", "HETATM", "ENDMDL", "CONECT", "MASTER", "END",
    "ENDMOL",
)
# MTRIX, ORIGX and SCALE pdb sections can be followed be a random
# number so they need to be looked at separately
PDB_SHORT_SECTION = ("MTRIX", "ORIGX", "SCALE")


class FormatError(Exception):
    """
    Exception raised when the file format is wrong.
    """
    pass


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


def valid_pdb_line(line):
    """
    Raise a FormatError if the given line does not start like a valid PDB line
    """
    if not (line[0:6].strip() in PDB_SECTIONS or line.strip() == ""):
        # MTRIX, ORIGX and SCALE pdb sections can be followed be a random
        # number so they will trigger the previous test
        if not line[0:4] in PDB_SHORT_SECTION:
            raise FormatError('PDB line should not start with "{0}"'
                              .format(line[0:6]))


def parse_selection(selection):
    """
    Read the atom selection given in argument

    The atom selection is formatted as follow:

    ::
        POPC:PO4 DUPC:PO4 CHOL:ROH

    Each string separated by a space represents one atom type, before the
    column is the residue name, after it is the atom name. The residue name can
    be omitted, then the column is omitted too:

    ::
        PO4 CHOL:ROH

    The final output is a list of tuples. The each tuple represents an atom
    type, the first element of the tuple is the residue name, the second
    element is the atom name. When the residue name is omitted then the tuple
    counts only one element: the atom name.

    Because the function is called from argparse the split on spaces is already
    done before entering the function, so here we only deal with one single
    atom type.
    """
    return tuple(selection.split(':'))


def stop_at_empty_line(iterator):
    """
    Yield all item of an iterator but stop when the item is an empty line.

    An empty line is a string which is empty when stripped.
    """
    for line in iterator:
        if line.strip() == "":
            return
        yield line


def except_last(iterator):
    """
    Yield all elements of an iterator but the last one.

    :Parameters:
        - iterator: the iterator on which to iterate
    """
    previous = next(iterator)
    for line in iterator:
        yield previous
        previous = line


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
    # "lines" might be a list and not a proper iterator
    lines = iter(lines)
    # The two first lines are a header
    next(lines)
    next(lines)
    # Loop over the lines, stop before an empty line and ignore the last
    # non empty line since it describes the box
    for line in except_last(stop_at_empty_line(lines)):
        try:
            atom = dict(((key, convert(line[begin:end].strip()))
                        for key, ((begin, end), convert)
                        in GRO_FIELDS.items()))
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
        valid_pdb_line(line)
        if line[0:6] == "ATOM  ":
            try:
                atom = dict(((key, convert(line[begin:end].strip()))
                            for key, ((begin, end), convert)
                            in PDB_FIELDS.items()))
            except ValueError:
                raise FormatError
            yield atom


def is_selected(atom, selection):
    """
    Return True is the atom fit the selection criteria.
    """
    for atom_type in selection:
        if len(atom_type) == 1 and atom["atom_name"] == atom_type[0]:
            return True
        if (atom["resname"] == atom_type[0] and
                atom["atom_name"] == atom_type[1]):
            return True
    return False


def select_atoms(atoms, selection):
    """
    Select only the atoms with the given atom name and residue name.

    This function create a generator that yield the dictionary for the atoms
    with the given atom name and residue name.

    :Parameters:
        - atoms: an iterator over the atom dictionaries
        - selection: an atom selection as described in ``parse_selection``
    """
    for atom in atoms:
        if is_selected(atom, selection):
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


def split_get_res(atoms, average, axis, selection):
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
        if is_selected(atom, selection) and not keep_res:
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
    for group_name, atomids in groups.items():
        print("[ {0} ]".format(group_name))
        group_str = " ".join([str(i) for i in atomids])
        print("\n".join(textwrap.wrap(group_str, 80)))


def split_leaflets(infile, axis, selection, res=False, input_format="gro"):
    """
    Split bilayer leaflets from a gromacs gro file along the given axis.
    """
    axis = axis.lower()
    if input_format == "gro":
        atoms = list(read_gro(infile))
    else:
        atoms = list(read_pdb(infile))
    selected = list(select_atoms(atoms, selection))
    coordinates = axis_coordinates(selected, axis)
    average = mean(coordinates)
    if res:
        groups = split_get_res(atoms, average, axis, selection)
    else:
        groups = split(selected, average, axis)
    write_ndx(groups)
    return groups


def get_options(argv):
    """
    Read the command line arguments.
    """
    usage = ("%(prog)s [options] < input > output.ndx\n"
             "       %(prog)s [options] input > output.ndx")
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
    keep_options.add_argument("--keep_residue", "-r", action="store_true",
                              dest="keep_residue", default=False,
                              help="Keep the whole residues.")
    keep_options.add_argument("--keep_atom", "-k", action="store_false",
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
    read_lines = []
    # Empty lines are not informative, let's go to the first not empty line
    line = infile.readline()
    while len(line) >= 1 and line.strip() == "":
        read_lines.append(line)
        line = infile.readline()
    read_lines.append(line)

    # If the file is empty it is not worth continuing
    if not line:
        return 'empty', []

    # Check if the line could be from a PDB file, if not we probably are
    # reading a gro file
    try:
        valid_pdb_line(line)
    except FormatError:
        input_format = "gro"
    else:
        input_format = "pdb"
    mod_infile = itertools.chain(read_lines, infile)
    return input_format, mod_infile


def reformat_selection(selection):
    """
    Generate a human readable string from an atom selection critera list.
    """
    return " ".join([":".join(criterion) for criterion in selection])


def main():
    """
    Run everything from the command line.
    """
    args = get_options(sys.argv[1:])
    if args.input is None:
        infile = sys.stdin
    else:
        try:
            infile = open(args.input)
        except IOError as error:
            print("Error while oppening file {0}".format(error.filename),
                  file=sys.stderr)
            return 1
    with infile:
        # Guess the format
        if args.format == "auto":
            input_format, mod_infile = guess_format(infile)
        else:
            input_format = args.format
            mod_infile = infile
        # Complain if the file is known to be empty
        if input_format == 'empty':
            print("The file is empty!", file=sys.stderr)
            return 1
        # Do the work
        try:
            groups = split_leaflets(mod_infile, args.axis, args.atom,
                                    args.keep_residue, input_format)
        # Complain if the format is wrong
        except FormatError:
            if (args.format == "auto"):
                print("Error while reading the input. Are you sure your file "
                      "is in the pdb or gro format?", file=sys.stderr)
            else:
                print(("Error while reading the input. Are you sure your file "
                       "is in the {0} format?").format(args.format),
                      file=sys.stderr)
                return 1
        # Complain if the reference atom is absent
        except ZeroDivisionError:
            print(("There seems to be no atom corresponding to your reference "
                   "selection. Are you sure this selection is present in your "
                   "structure: '{0}'?")
                  .format(reformat_selection(args.atom)), file=sys.stderr)
            return 1
    # Display the number of atoms per group
    for group_name, atomids in groups.items():
        print("{0}: {1} atoms".format(group_name, len(atomids)),
              file=sys.stderr)
    if len(groups["upper_leaflet"]) == len(groups["lower_leaflet"]):
        print("The membrane is symmetric.", file=sys.stderr)

    return 0

if __name__ == "__main__":
    sys.exit(main())

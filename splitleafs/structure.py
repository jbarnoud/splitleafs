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
Parse structure files in gro and pdf format.
"""

from __future__ import print_function, division, with_statement
import itertools

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



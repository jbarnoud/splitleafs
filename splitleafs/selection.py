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
Select atoms.
"""

from __future__ import print_function, division, with_statement

__author__ = "Jonathan Barnoud"


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


def reformat_selection(selection):
    """
    Generate a human readable string from an atom selection critera list.
    """
    return " ".join([":".join(criterion) for criterion in selection])



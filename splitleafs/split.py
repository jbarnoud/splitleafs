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
from .selection import is_selected, select_atoms

__author__ = "Jonathan Barnoud"


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


def split_leaflets(infile, axis, selection, file_reader, res=False):
    """
    Split bilayer leaflets from a gromacs gro file along the given axis.

    :Parameters:
        - infile: the input file describing the structure as a iterator over
                  the structure file (gro or pdb format)
        - axis: the dimension normal to the membrane plane (x, y, or z)
        - selection: an atom selection list as outputed by ``parse_selection'',
                     the selection is used to get the reference atoms
        - file_reader: a callback to the function that will read the input
                       (typically read_gro or read_pdb)
        - res: a boolean, True if you want to keep whole residues in the output,
               False by default

    :Return:
        - a dictionary, they keys are "upper_leaflet" and "lower_leaflet", the
          values are lists of atom indices in each leaflet. The indices start
          at 1!
    """
    axis = axis.lower()
    atoms = list(file_reader(infile))
    selected = list(select_atoms(atoms, selection))
    coordinates = axis_coordinates(selected, axis)
    average = mean(coordinates)
    if res:
        groups = split_get_res(atoms, average, axis, selection)
    else:
        groups = split(selected, average, axis)
    return groups



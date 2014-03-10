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
import sys
from functools import partial

from .structure import read, FormatError
from .user import get_options
from .split import split_leaflets
from .output import write_ndx

__author__ = "Jonathan Barnoud"


def main():
    """
    Run everything from the command line.
    """
    args = get_options(sys.argv[1:])
    if args.input is None:
        print("Read input from the standard input.", file=sys.stderr)
        infile = sys.stdin
    else:
        try:
            infile = open(args.input)
        except IOError as error:
            print("Error while oppening file {0}".format(error.filename),
                  file=sys.stderr)
            return 1
    with infile:
        file_reader = partial(read, file_format=args.format)
        # Do the work
        try:
            groups = split_leaflets(infile, args.axis, args.atom,
                                    file_reader, args.keep_residue)
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
        else:
            write_ndx(groups)

    # Display the number of atoms per group
    for group_name, atomids in groups.items():
        print("{0}: {1} atoms".format(group_name, len(atomids)),
              file=sys.stderr)
    if len(groups["upper_leaflet"]) == len(groups["lower_leaflet"]):
        print("The membrane is symmetric.", file=sys.stderr)

    return 0

if __name__ == "__main__":
    sys.exit(main())

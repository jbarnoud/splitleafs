Splitleaf
=========

This tools aims to split leaflets from a bilayer. It reads a GROMACS gro file
and print an index file with one group per leaflet.

Usage::

    splitleafs.py --axis z --atom P1 < input.gro > output.ndx

Options:

--axis or -d: The axis normal to the membrane, z by default
--atom or -a: The atom name of the reference, P1 by default

Splitleafs
==========

.. image:: https://travis-ci.org/jbarnoud/splitleafs.png?branch=master
   :alt: Travis-CI build status
   :target: https://travis-ci.org/jbarnoud/splitleafs

This tools aim to split leaflets from a bilayer. It reads a GROMACS gro file or
a PDB file and print an index file with one group per leaflet.

Usage
-----
::

    splitleafs.py [options] < input > output.ndx
    splitleafs.py [options] input > output.ndx


If an input file if given as an argument then this file will be read; else the
program will read the standard input. The input file can be in the gro or PDB
format. The format is guessed from the content of the file Use the
``--format`` if you want to specify explicitly the file format.

The program writes a gromacs index file (.ndx file) on the standard output. The
index file will contain two groups called ``lower_leaflet`` and
``upper_leaflet`` respectively. The groups contain the index of the reference
atoms in each leaflets (if the ``--keep_atom`` option is used, default) or
the index of all the atoms from reference atom residues (if the
``--keep_residue`` is used).

Splitleafs has been tested with python 2.6, 2.7, 3.2 and pypy.

Options
~~~~~~~

--axis or -d:
    The axis normal to the membrane, ``z`` by default.
--atom or -a:
    The atoms of reference. Several Several atoms can be given. The format for
    one reference atom is ``residue_name:atom_name``, the residue name can be
    omitted then the atom is mentioned by its atom name only.
--keep_residue or -r:
    Write the whole residues in the index file.
--keep_atom or -k:
    Write only the atoms of reference in the index file (default).
--format or -f:
    The input file format. Gromacs gro files and PDB files are supported.
    If ``auto`` is selected, then the file format is guessed from the input
    content. ``auto`` by default.

Examples
~~~~~~~~
::
    
    splitleafs.py membrane.pdb > leaflets.ndx

This command will read the ``membrane.pdb`` file, it will guess it is a PDB file
and will write the ``leaflets.ndx`` file with the index of the ``P1`` atoms
split in two leaflets. Here we assume that ``membrane.pdb`` describe a
phospholipid and that ``P1`` is the name of the phosphorus.

If your system is different you can change the reference atom using the
``--atom`` option. Let's assume we have a phospholipid bilayer described within
the `Martini coarse-grained force field <http://md.chem.rug.nl/cgmartini/>`_.
Then the phospate group is represented by a bead called ``PO4``. If we want to
use this bead as a reference for our splitting we can use the command ::

    splitleafs.py --atom PO4 membrane_martini.gro > leaflets.ndx

Notice that, in this last example, we give a gro file to the program which will
guess it is not a PDB file and read it with the right file format.

If the membrane in the previous example is made of POPC lipids only, then the
command is equivalent to ::

    splitleafs.py --atom POPC:PO4 membrane_martini.gro > leaflets.ndx

The ``--atom`` option can take several arguments. This can be useful in case of
an heterogeneous membrane ::

    splitleafs.py --atom POPC:PO4 CHOL:ROH membrane_martini.gro > leaflets.ndx

You may need not only the phosphorus but also the whole lipid in the output
index file. Then use the ``--keep_residue`` option ::

    splitleafs.py --atom PO4 membrane_martini.gro --keep_residue > leaflets.ndx

In the current version, only the residue that contains the reference atom is
selected. This means that if your molecules are composed of more than one
residue, the others will be ignored. Most lipid representation use only one
residue per lipid so this issue should not be a concern.

You might want to know how many atoms were selected in each leaflet. This is
convenient to quickly check the result of the program or to know if a membrane
is symmetric. These informations are written on the standard error stream so you
can see them even if you redirect the output of the program in a file ::

    lower_leaflet: 144 atoms
    upper_leaflet: 144 atoms
    The membrane is symmetric.


License
-------

This program is free software: you can redistribute it and/or modify  
it under the terms of the GNU General Public License as published by   
the Free Software Foundation, either version 3 of the License, or      
(at your option) any later version.                                    
                                                                      
This program is distributed in the hope that it will be useful,        
but WITHOUT ANY WARRANTY; without even the implied warranty of         
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
GNU General Public License for more details.                           
                                                                          
A copy of the GNU General Public License is available at
http://www.gnu.org/licenses/gpl-3.0.html.


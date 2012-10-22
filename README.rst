Splitleaf
=========

This tools aims to split leaflets from a bilayer. It reads a GROMACS gro file
and print an index file with one group per leaflet.

Usage
-----
::

    splitleafs.py --axis z --atom P1 < input.gro > output.ndx

Options:

--axis or -d: The axis normal to the membrane, z by default
--atom or -a: The atom name of the reference, P1 by default
--keep_residue or -r: Write the whole residues in the index file
--keep_atom or -k: Write only the atoms of reference in the index file (default)

Licence
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


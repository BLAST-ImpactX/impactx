#!/usr/bin/env python3

import sys, os
import numpy as np
import synergia

lattice_madx =     fodo_madx = """
beam, particle=proton,energy=pmass+0.00250;
dedge30: dipedge,e1:= 0,h:= 1.216203238,hgap:= 0.029,fint:= 0.5;

lattice: sequence, l=0.0, refer=entry;
dedge30, at=0.0;
endsequence;
"""

reader = synergia.lattice.MadX_reader()
reader.parse(lattice_madx)
lattice = reader.get_lattice('lattice')

print(lattice)

refpart = lattice.get_reference_particle()
print('energy: ', refpart.get_total_energy())

map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)
print(map)

print('R43 element: ', map[3,2])

# value from running madx script print_dipedge_map.madx in file mymap.out
madx_43 = 0.04291315480354536
print('madx R43: ', madx_43)

#!/usr/bin/env python
import sys
import numpy as np
import synergia as syn

# print the third order map for a lattice

# the coordinate names
cnm = ["x", "xp", "y", "yp", "cdt", "dpop" ]

# given a lattice, print it's 3rd order maps

def print_o3(lattice, ofile=sys.stdout):

    # calculate one turn map at given order (3 here)
    mapping = syn.simulation.Lattice_simulator.get_one_turn_map_o3(lattice)

    #  iterate through the components and fields
    for comp in range(6):
        print(f'\nCoordinate {cnm[comp]} = ', file=ofile)
        trigon = mapping.component(comp)

        for pwr in range(trigon.power()+1):
           for idx in range(trigon.count(pwr)): 
               idx_to_exp = "syn.foundation.Trigon_index_to_exp_o{}(idx)".format(pwr)
               coeff = trigon.get_term(pwr, idx)
               if abs(coeff) < 1.0e-15:
                   continue

               print(f'{coeff} ', end='', file=ofile)

               # get the exponents of the polynomial term
               poly_exps = eval(idx_to_exp)
               
               for c,p in enumerate(poly_exps):
                   if p != 0:
                       print(f'{cnm[c]}^{p} ', end='', file=ofile)

               print(file=ofile)

usage_txt = """
print_map_o3.py <lattice-file> <seq-name>

prints the 3rd order maps for the sequence <seq-name> within the
lattice file given by <lattice-file>. There must be a BEAM statement
within the file to specify the particle mass and energy.
"""

def main(argv):
    if len(argv) < 3:
        print(usage_txt)
        sys.exit(10)

    reader = syn.lattice.MadX_reader()
    lattice = reader.get_lattice(argv[2], argv[1])
    syn.simulation.Lattice_simulator.tune_circular_lattice(lattice)

    print_o3(lattice)

if __name__ == "__main__":
    main(sys.argv)

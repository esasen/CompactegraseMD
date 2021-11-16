#!/bin/env python3

import os,sys
from typing import List, Tuple

import gen_DNA

try:
    import numpy as np
except ModuleNotFoundError:
    print('numpy not installed. Please install numpy: pip install numpy')
    sys.exit(1)

########################################################################
########################################################################
########################################################################

class LmpSetup:

    dimensions  = 3
    units       = 'lj'
    boundary    = 'ppp'
    temp        = 1
    gamma       = 1
    atom_style  = 'angle'

    seed: int

    def __init__(self,simbox,seed=-1):
        if seed == -1:
            self.gen_seed()
            print(self.seed)

    ######################################
    # seed functions
    def gen_seed(self):
        self.seed = np.random.randint(0,1000000000)


########################################################################
########################################################################
########################################################################

class AtomType:
    name:       str
    id:         int

    radius:     float
    mass:       float

    def __init__(self,name,mass,radius,id=-1):
        self.name       = name
        self.mass       = mass
        self.radius     = radius
        self.id         = id

class BondTypes:
    name:       str
    id:         int

    bond_style: str
    bond_coeff: List

    def __init(self,name,bond_style,bond_coeff):
        self.bond_style = bond_style
        self.bond_coeff = bond_coeff

class AngleTypes:
    name:       str
    id:         int

    bond_style: str
    bond_coeff: List

    def __init(self,name,bond_style,bond_coeff):
        self.bond_style = bond_style
        self.bond_coeff = bond_coeff




# class MoleculeType:
#
#     contained_atoms:    List[atom_type]
#     fixed_atom:         bool
#     bond_pairs:         List[atom_type]
#
#     def __init__(self,atoms):
#         self.atoms = atoms


# class Bonds:
#
#
#     def __init__(self,id1,id2):








########################################################################
########################################################################
########################################################################

if __name__ == "__main__":
    
    # if len(sys.argv) < 2:
    #     print("usage: %s num_bp L"%sys.argv[0])
    #     sys.exit()
    
    # nbp         = int(sys.argv[1])
    # disc_len    = float(sys.argv[2])
    # L           = float(sys.argv[3])

    L = 200
    simbox = np.array([[0, L], [0, L], [0, L]])

    lmp = LmpSetup(simbox=simbox)
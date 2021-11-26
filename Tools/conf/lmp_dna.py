#!/bin/env python3

import os,sys
from typing import List, Tuple

import gen_DNA
import lmp_setup as lmps

try:
    import numpy as np
except ModuleNotFoundError:
    print('numpy not installed. Please install numpy: pip install numpy')
    sys.exit(1)

########################################################################
########################################################################
########################################################################

class DNA(lmps.molecule_type):







class AtomType:
    name:       str
    mol_name:   str

    id:         int
    molid:      int


    radius:     float
    mass:       float
    pos :       np.ndarray

    def __init__(self,name,mass,radius,mol_name=None):
        self.name       = name
        self.mass       = mass
        self.radius     = radius
        self.mol_name   = mol_name


class MoleculeType:

    contained_atoms: List[atom_type]
    fixed_atom:      bool

    def __init__(self,atoms):
        self.atoms = atoms


class Bonds:









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

    lmp = lmp_setup(simbox=simbox)
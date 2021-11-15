#!/bin/env python3

import sys
try:
    import numpy as np
except ModuleNotFoundError:
    print('numpy not installed. Please install numpy')
    sys.exit(1)

try:
    from numba import jit
    USE_NUMBA = True
except ModuleNotFoundError:
    print('Warning: numba not installed. Consider installing numba to speed up the calculation')
    USE_NUMBA = False

import

########################################################################
########################################################################
########################################################################

@jit(nopython=True)
def ev_violation_pair(r1,r2,ev_dist):
    """ checks if pair of atoms violates excluded volume set by ev_dist """




if __name__ == "__main__":

    box = np.array([[-2.,10],[-2.,10],[-2,10]])
    # box = np.array([[-2.,10],[-2,10]])


    r = np.array([[5,80,-9],[-456,46.57,-2.001]])
    print(r)
    r = place_in_box(box,r)

    print(r)




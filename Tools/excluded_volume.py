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

########################################################################
########################################################################
########################################################################

@jit(nopython=True)
def ev_violation_pair(r1,r2,ev_dist: float):
    if np.linalg.norm(r2-r1) <= ev_dist:
        return False
    return True

@jit(nopython=True)
def ev_violation_pair_relsize(r1,r2,size1: float,size2: float):
    ev_dist = 0.5*(size1+size2)
    return ev_violation_pair(r1, r2, ev_dist)

@jit(nopython=True)
def ev_violation(R,ev_dist,excluded_neighbors=0: int):
    """ checks for excluded volume violations within the given group of atoms """


    ev_dist = 0.5*(size1+size2)
    return ev_violation_pair(r1, r2, ev_dist)


def place_in_box_single(box,r):
    """ places a single position within the bounds of the specified box. The box dimension is 3x2 """
    if np.shape(box) != (3,2):
        raise ValueError(f'Invalid dimension of box. Needs to be (3x2), {np.shape(box)} given')

    # use numba is numba is loaded
    if USE_NUMBA:
        return place_in_box_single_numba(box,r)

    rnew = np.zeros(np.shape(r))
    rnew[0] = (r[0]-box[0,0])%(box[0,1]-box[0,0]) + box[0,0]
    rnew[1] = (r[1]-box[1,0])%(box[1,1]-box[1,0]) + box[1,0]
    rnew[2] = (r[2]-box[2,0])%(box[2,1]-box[2,0]) + box[2,0]
    return rnew

if USE_NUMBA:
    @jit(nopython=True)
    def place_in_box_single_numba(box,r):
        rnew = np.zeros(np.shape(r))
        rnew[0] = (r[0]-box[0,0])%(box[0,1]-box[0,0]) + box[0,0]
        rnew[1] = (r[1]-box[1,0])%(box[1,1]-box[1,0]) + box[1,0]
        rnew[2] = (r[2]-box[2,0])%(box[2,1]-box[2,0]) + box[2,0]
        return rnew

########################################################################
# impose periodic boundary for array of atoms

def place_in_box(box,R):
    """ places an array of positions within the bounds of the specified box. The box dimension is 3x2 """
    if np.shape(box) != (3,2):
        raise ValueError(f'Invalid dimension of box. Needs to be (3x2), {np.shape(box)} given')

    # use numba is numba is loaded
    if USE_NUMBA:
        return place_in_box_numba(box,R)

    Rnew = np.zeros(np.shape(R))
    span = box[:,1] - box[:,0]
    for i in range(len(R)):
        Rnew[i,0] = (R[i,0]-box[0,0])%span[0] + box[0,0]
        Rnew[i,1] = (R[i,1]-box[1,0])%span[1] + box[1,0]
        Rnew[i,2] = (R[i,2]-box[2,0])%span[2] + box[2,0]
    return Rnew

if USE_NUMBA:
    @jit(nopython=True)
    def place_in_box_numba(box,R):
        Rnew = np.zeros(np.shape(R))
        span = box[:,1] - box[:,0]
        for i in range(len(R)):
            Rnew[i,0] = (R[i,0]-box[0,0])%span[0] + box[0,0]
            Rnew[i,1] = (R[i,1]-box[1,0])%span[1] + box[1,0]
            Rnew[i,2] = (R[i,2]-box[2,0])%span[2] + box[2,0]
        return Rnew

def valid_box_dimension(box):
    if np.shape(box) != (3,2):
        raise ValueError(f'Invalid dimension of box. Needs to be (3x2), {np.shape(box)} given')

def valid_box(box):
    if np.shape(box) != (3,2):
        raise ValueError(f'Invalid dimension of box. Needs to be (3x2), {np.shape(box)} given')
    for i in range(3):
        if box[i,1] <= box[i,0]:
            raise ValueError(f'Lower bound of periodic box larger than upper bound')


if __name__ == "__main__":

    box = np.array([[-2.,10],[-2.,10],[-2,10]])
    # box = np.array([[-2.,10],[-2,10]])


    r = np.array([[5,80,-9],[-456,46.57,-2.001]])
    print(r)
    r = place_in_box(box,r)

    print(r)




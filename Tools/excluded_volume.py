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
def ev_violation_pair(r1,r2,ev_dist):
    """ checks if pair of atoms violates excluded volume set by ev_dist """
    if np.linalg.norm(r2-r1) < ev_dist:
        return True
    return False

@jit(nopython=True)
def ev_violation_pair_relsize(r1,r2,size1: float,size2: float):
    """ checks if pair of atoms violates excluded volume set by the respective size of the atoms """
    ev_dist = 0.5*(size1+size2)
    return ev_violation_pair(r1, r2, ev_dist)

@jit(nopython=True)
def ev_violation(R,ev_dist,excluded_neighbors=0):
    """ checks for excluded volume violations within the given group of atoms """
    N = len(R)
    for i in range(excluded_neighbors+1,N):
        for j in range(0,i-excluded_neighbors):
            if np.linalg.norm(R[j] - R[i]) < ev_dist:
                return True
    return False

@jit(nopython=True)
def ev_violation_single(R,r,ev_dist):
    """ checks if the atom specified by r violates excluded volumes with the group of atoms R """
    for i in range(len(R)):
        if np.linalg.norm(r - R[i]) < ev_dist:
            return True
    return False

########################################################################
# including periodic boundary conditions

@jit(nopython=True)
def ev_violation_single_periodic_box(R,r,ev_dist,periodic_box):
    """ checks if the atom specified by r violates excluded volumes with the group of atoms R including periodic boundary"""
    N = len(R)
    for i in range(N):
        if ev_violation_pair_in_periodic_box(r,R[i],ev_dist,periodic_box):
            return True
    return False


@jit(nopython=True)
def ev_violation_pair_in_periodic_box(r1,r2,ev_dist,periodic_box):
    """
     checks if the two atoms or their periodic copies violate the excluded volume set by ev_dist
     assumes the tho positions to be within the box
    """
    dists = np.zeros(3)
    for i in range(3):
        dists[i] = closest_dist_1d_periodic(r1[i], r2[i], periodic_box[i])
    if np.linalg.norm(dists) < ev_dist:
        return True
    return False

@jit(nopython=True)
def closest_dist_1d_periodic(a,b,box_bound):
    """ calculates the closest distance of two points in 1d across the periodic boundary """
    dx = np.abs(a - b)
    box_range = box_bound[1] - box_bound[0]
    if a > box_bound[0] + box_range*0.5:
        dx_alt = np.abs(a - box_range - b)
    else:
        dx_alt = np.abs(a + box_range - b)
    if dx_alt < dx:
        return dx_alt
    return dx



if __name__ == "__main__":

    box = np.array([[-2.,10],[-2.,10],[-2,10]])
    # box = np.array([[-2.,10],[-2,10]])


    r = np.array([[5,80,-9],[-456,46.57,-2.001]])
    print(r)
    r = place_in_box(box,r)

    print(r)




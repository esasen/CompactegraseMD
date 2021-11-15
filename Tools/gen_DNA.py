#!/bin/env python3

import os,glob,sys
import random
import SO3Methods as so3
import periodic_boundary

try:
    import numpy as np
except ModuleNotFoundError:
    print('numpy not installed. Please install numpy: pip install numpy')
    sys.exit(1)

try:
    from numba import jit
except ModuleNotFoundError:
    print('numba not installed. Please install numpy: pip install numba')
    sys.exit(1)

    
########################################################################
########################################################################
########################################################################

def gen_DNA_conf(nbp,disc_len,periodic_box=None,lb=40,excluded_volume=0,first_pos=None,first_triad=None):
    """ generates a DNA configuration. periodic boundaries will be imposed if the periodic_box argument is passed.
        Excluded columes between atoms will be imposed based on an excluded volume radius set by the argument.
    """

    # check if periodic boundary box dimensions are correct
    if periodic_box is not None:
        periodic_boundary.valid_box(periodic_box)

    if first_pos is None:
        if periodic_box is None:
            first_pos = np.zeros(3)

    if first_triad is None:
        theta = np.random.random(3) * np.pi
        first_triad = so3.get_rot_mat(theta)

    pos = gen_conf(nbp,disc_len,lb,first_pos,first_triad)

    if periodic_box is not None:
        pos = periodic_boundary.place_in_box(periodic_box,pos)




@jit(nopython=True)
def gen_conf(nbp,disc_len,lb,first_pos,first_triad):
    pos    = np.zeros((nbp,3))
    T      = np.copy(first_triad)
    sigma = np.sqrt(disc_len/lb)

    pos[0] = first_pos
    pos[1] = first_pos + T[:,2]*disc_len
    for i in range(2,nbp):
        theta = np.random.normal(loc=0.0, scale=sigma, size=3)
        R = so3.get_rot_mat(theta)
        T = np.dot(T,R)
        pos[i] = pos[i-1] + T[:, 2] * disc_len
    return pos











if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("usage: %s num_bp L"%sys.argv[0])
        sys.exit()
    
    nbp         = int(sys.argv[1])
    disc_len    = float(sys.argv[2])
    L           = float(sys.argv[3])

    # theta = np.array([0,np.pi/2,0])
    # theta = np.array([0.,0,0])
    # R = so3.get_rot_mat(theta)
    # print(R)
    # print(R[:,2])

    gen_DNA_conf(nbp, disc_len, periodic_box=None, lb=40, excluded_volume=0, first_pos=None, first_triad=None)

    
    # for i in range(nbp):
    #     if i%10000 == 0:
    #         print(i)
    #
    #     theta = np.random.random(3)*0.99*np.pi
    #     # ~ print(theta)
    #     R = so3.get_rot_mat(theta)
        

#!/bin/env python3

import sys
import SO3Methods as so3
import periodic_boundary as bp
import excluded_volume as ev

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

from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

    
########################################################################
########################################################################
########################################################################

def gen_DNA_conf(nbp: int,disc_len: float,periodic_box=None,lb=40.0,ev_size=0.0,first_pos=None,first_triad=None,rise=0.34):
    """ generates a DNA configuration. periodic boundaries will be imposed if the periodic_box argument is passed.
        Excluded columes between atoms will be imposed based on an excluded volume radius set by the argument.
    """

    num_segs = int(np.ceil(nbp*rise/disc_len))

    # check if periodic boundary box dimensions are correct
    if periodic_box is not None:
        bp.valid_box(periodic_box)

    # set first pos
    if first_pos is None:
        if periodic_box is None:
            first_pos = np.zeros(3)
        else:
            first_pos = gen_random_point_in_box(periodic_box)

    # set first triad
    if first_triad is None:
        theta = np.random.random(3) * np.pi
        first_triad = so3.get_rot_mat(theta)

    if ev_size is None or ev_size == 0:
        # if no excluded volume, no trials are required
        pos = gen_conf(num_segs, disc_len, lb, first_pos, first_triad)
        if periodic_box is not None:
            pos = bp.place_in_box(periodic_box,pos)
    else:
        # with excluded volume the somewhat more expensive generation function is required
        excluded_neighbors = int(np.ceil(ev_size / disc_len))
        pos = gen_ev_conf(num_segs,
                          disc_len,
                          lb,
                          first_pos,
                          first_triad,
                          ev_size,
                          excluded_neighbors=excluded_neighbors,
                          periodic_box=periodic_box)
    return pos


# @jit(nopython=True)
def gen_random_point_in_box(periodic_box):
    return periodic_box[:,0] + np.random.rand(3)*(periodic_box[:,1]-periodic_box[:,0])

@jit(nopython=True)
def gen_conf(num_segs,disc_len,lb,first_pos,first_triad):
    """
        Generate configuration without taking excluded volumes into account
    """
    pos    = np.zeros((num_segs,3))
    T      = np.copy(first_triad)
    sigma = np.sqrt(disc_len/lb)

    pos[0] = first_pos
    pos[1] = first_pos + T[:,2]*disc_len
    for i in range(2,num_segs):
        theta = np.random.normal(loc=0.0, scale=sigma, size=3)
        R = so3.get_rot_mat(theta)
        T = np.dot(T,R)
        pos[i] = pos[i-1] + T[:, 2] * disc_len
    return pos

# @jit(nopython=True)
def gen_ev_conf(num_segs,disc_len,lb,first_pos,first_triad,ev_size,excluded_neighbors=0,periodic_box=None,shift_back = 20,max_trials_per_step=1000):
    """
        Generate configuration considering exlcluded volumes and the periodicity of the box
    """
    pos         = np.zeros((num_segs,3))
    triads      = np.zeros((num_segs-1,3,3))
    # T           = np.copy(first_triad)5
    triads[0]   = first_triad
    sigma = np.sqrt(disc_len/lb)

    pos[0] = first_pos
    pos[1] = first_pos + first_triad[:,2]*disc_len
    i = 2

    max_trials = max_trials_per_step*num_segs
    trials     = 0

    while i < num_segs:
        trials += 1
        if trials >= max_trials:
            raise Exception(f"Could not generate configuration. Max trial steps ({max_trials}) reached")

        theta = np.random.normal(loc=0.0, scale=sigma, size=3)
        # print(theta)
        R = so3.get_rot_mat(theta)
        triads[i-1] = np.dot(triads[i-2],R)
        pos[i] = pos[i-1] + triads[i-1,:, 2] * disc_len

        # if periodic box is not active
        if periodic_box is None:
            if (i - excluded_neighbors) > 0:
                if ev.ev_violation_single(pos[:i - excluded_neighbors], pos[i], ev_size):
                    print(f'reverting at step {i}')
                    # violated
                    i = i - shift_back
                    if i < 2:
                        i = 2
                    continue

        # if periodic box is active
        else:
            pos[i] = bp.place_in_box_single(periodic_box, pos[i])
            if (i - excluded_neighbors) > 0:
                if ev.ev_violation_single_periodic_box(pos[:i - excluded_neighbors], pos[i], ev_size, periodic_box):
                    # violated
                    i = i - shift_back
                    if i < 2:
                        i = 2
                    continue
        i+=1
    return pos

def plot_DNA_conf(conf):
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    ax.scatter(conf[:,0], conf[:,1], conf[:,2], zdir='z', s=50, c='black', depthshade=True)
    ax.plot(conf[:, 0], conf[:, 1], conf[:, 2], zdir='z', c='black',lw=3)
    plt.show()



if __name__ == "__main__":
    
    if len(sys.argv) < 3:
        print("usage: %s num_bp disc_len L"%sys.argv[0])
        sys.exit()
    
    nbp         = int(sys.argv[1])
    disc_len    = float(sys.argv[2])
    L           = float(sys.argv[3])

    ev_size = disc_len
    ev_size = 0
    periodic_box = np.array([[0,L],[0,L],[0,L]])
    # periodic_box = None


    conf = gen_DNA_conf(nbp, disc_len, periodic_box=periodic_box, lb=40, ev_size=ev_size, first_pos=None, first_triad=None)

    # check ev
    for i in range(len(conf)):
        for j in range(i-2):
            if periodic_box is None:
                if np.linalg.norm(conf[i]-conf[j]) < ev_size:
                    print('EV VIOLATION!')
            else:

                if ev.ev_violation_pair_in_periodic_box(conf[i],conf[j],ev_size,periodic_box):
                    print('EV VIOLATION!')

    plot_DNA_conf(conf)

    # uconf = bp.unwrap_dna(conf,periodic_box,disc_len=None)
    # plot_DNA_conf(uconf)

    check_unwrap = False
    if check_unwrap:
        M = 20
        confs = np.zeros([M,nbp,3])
        for m in range(M):
            print(m)
            confs[m] = gen_DNA_conf(nbp, disc_len, periodic_box=periodic_box, lb=40, ev_size=ev_size, first_pos=None, first_triad=None)

        uconfs = bp.unwrap_dna(confs, periodic_box, disc_len=None)
        for m,uconf in enumerate(uconfs):
            plot_DNA_conf(confs[m])
            plot_DNA_conf(uconfs[m])
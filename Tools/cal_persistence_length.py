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

import gen_DNA

########################################################################
########################################################################
########################################################################

def cal_persistence_length(conf,m_max=1,disc_len=None):
    """ calculates the persistence length for the given set of configurations """
    if len(np.shape(conf)) == 2:
        conf = np.array([conf])
    if disc_len is None:
        disc_len = cal_disc_len(conf)

    m_max = 1

    tans =  __get_tangents(conf,normalized=True)
    if m_max == 1:
        lb = __cal_perslen_single(tans,disc_len)
    else:
        lb = __cal_perslen_multiple(tans, m_max, disc_len)

    print(lb)


@jit(nopython=True)
def __cal_perslen_single(tans,disc_len):
    tancor = 0.
    for s in range(len(tans)):
        for i in range(len(tans[0])-1):
            tancor += np.dot(tans[s,i],tans[s,i+1])

    tancor /= len(tans)*(len(tans[0])-1)
    return -disc_len/np.log(tancor)


@jit(nopython=True)
def __cal_perslen_multiple(tans, m_max, disc_len):
    tancor = np.zeros(m_max)
    numcor = np.zeros(m_max)
    for s in range(len(tans)):
        for i in range(len(tans[0])-1):
            for m in range(m_max):
                j = i+1+m
                if j >= len(tans[0]):
                    break
                tancor[m] += np.dot(tans[s,i],tans[s,i+1+m])
                numcor[m] += 1
    return -disc_len*np.arange(1,m_max+1) / np.log(tancor/numcor)


def cal_disc_len(conf):
    """ returns the mean discretization lenght. """
    if len(np.shape(conf)) == 2:
        conf = np.array([conf])
    return np.round(__cal_disc_len(conf), decimals=8)

@jit(nopython=True)
def __cal_disc_len(conf):
    dlen = 0.0
    for s in range(len(conf)):
        for i in range(len(conf[0])-1):
            dlen += np.linalg.norm(conf[s,i+1]-conf[s,i])
    return dlen/(len(conf)*(len(conf[0])-1))


def get_tangents(conf, normalized=False):
    """ returns tangents for given configuration. Return matrix will be of dim (m,n,3) """
    if len(np.shape(conf)) == 2:
        conf = np.array([conf])
    return __get_tangents(conf,normalized=normalized)

@jit(nopython=True)
def __get_tangents(conf,normalized=False):
    tans = np.zeros((len(conf), len(conf[0]), 3))
    for s in range(len(conf)):
        for i in range(len(conf[0])-1):
            tans[s,i] = conf[s, i + 1] - conf[s, i]
            if normalized:
                tans[s, i] /= np.linalg.norm(tans[s,i])
    return tans







if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: %s num_bp disc_len" % sys.argv[0])
        sys.exit()

    nbp = int(sys.argv[1])
    disc_len = float(sys.argv[2])

    ev_size = None

    # periodic_box = np.array([[0, L], [0, L], [0, L]])
    periodic_box = None

    conf = gen_DNA.gen_DNA_conf(nbp, disc_len, periodic_box=periodic_box, lb=40, ev_size=ev_size, first_pos=None,
                        first_triad=None)
    print('done')
    lb = cal_persistence_length(conf, m_max=10, disc_len=None)



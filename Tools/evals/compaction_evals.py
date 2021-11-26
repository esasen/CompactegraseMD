#!/bin/env python3

import sys
try:
    import numpy as np
except ModuleNotFoundError:
    print('numpy not installed. Please install numpy')
    sys.exit(1)

try:
    from numba import jit
except ModuleNotFoundError:
    print('numba not installed. Please install numpy: pip install numba')
    sys.exit(1)


########################################################################
########################################################################
########################################################################

# @jit(nopython=True)
def radius_of_gyration(conf):7
    com  = center_of_mass(conf)
    return np.sqrt(np.sum((conf-com)**2))

def center_of_mass(conf):
    return np.mean(conf,axis=0)




if __name__ == "__main__":

    if len(sys.argv) < 4:
        print("usage: %s num_bp disc_len" % sys.argv[0])
        sys.exit()

    nbp = int(sys.argv[1])


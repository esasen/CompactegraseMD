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


def gen_DNA(nbp,disc_len,periodic_box=None,lb=40,excluded_volume=0):





    












if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("usage: %s num_bp L"%sys.argv[0])
        sys.exit()
    
    nbp = int(sys.argv[1])
    L   = float(sys.argv[2])
    
    
    for i in range(nbp):
        if i%10000 == 0:
            print(i)
        
        theta = np.random.random(3)*0.99*np.pi
        # ~ print(theta)
        R = so3.get_rot_mat(theta)
        

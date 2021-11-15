import os,glob,sys
import random
import SO3Methods as so3

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
        

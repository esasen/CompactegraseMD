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
# impose periodic boundary for single atom

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



if __name__ == "__main__":

    box = np.array([[-2.,10],[-2.,10],[-2,10]])
    # box = np.array([[-2.,10],[-2,10]])


    r = np.array([[5,80,-9],[-456,46.57,-2.001]])
    print(r)
    r = place_in_box(box,r)

    print(r)




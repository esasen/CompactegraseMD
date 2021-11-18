#!/bin/env python3

import sys

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

import gen_DNA as gdna


    
########################################################################
########################################################################
########################################################################

def gen_conffile_with_DNA(outfn: str, num_atom_types: int, num_bond_types: int, num_angle_types: int,periodic_box: np.ndarray, dna_conf,dna_types, dna_bond_type=1,dna_angle_type=1):

    with open(outfn,'w') as f:

        # write header
        f.write('HEADER\n\n')

        # write num atoms
        f.write('   %d atoms\n'%len(dna_conf))
        f.write('   %d bonds\n'%(len(dna_conf)-1))
        f.write('   %d angles\n'%(len(dna_conf)-2))

        f.write('\n')
        f.write('   %d atom types\n' % num_atom_types)
        f.write('   %d bond types\n' % num_bond_types)
        f.write('   %d angle types\n' % num_angle_types)

        f.write('\n')
        f.write('# Define periodic box\n')
        f.write('   %.8f %.8f xlo xhi\n' % (periodic_box[0,0], periodic_box[0,1] ))
        f.write('   %.8f %.8f ylo yhi\n' % (periodic_box[1,0], periodic_box[1,1] ))
        f.write('   %.8f %.8f zlo zhi\n' % (periodic_box[2,0], periodic_box[2,1] ))

        f.write('\n')
        f.write('Atoms\n\n')
        for i in range(len(dna_conf)):
            type     = 1
            molecule = 1
            f.write('%d %d %d %.8f %.8f %.8f\n'%(i+1,type,molecule,dna_conf[i,0],dna_conf[i,1],dna_conf[i,2]))

        f.write('\n')
        f.write('Bonds\n\n')
        for i in range(len(dna_conf)-1):
            f.write('%d %d %d %d\n'%(i+1,dna_bond_type,i+1,i+2))

        f.write('\n')
        f.write('Angles\n\n')
        for i in range(len(dna_conf)-2):
            f.write('%d %d %d %d %d\n'%(i+1,dna_angle_type,i+1,i+2,i+3))


# def __get_dna_type_list(num_segs,dna_types):
#     if len(dna_types) == 1:
#         types =

#
# Test Conf
#
# 	0 atoms
#
# 	3 atom types
#
#
# 	-50	50 xlo xhi
# 	-50	50 ylo yhi
# 	-50	50 zlo zhi



if __name__ == "__main__":
    
    if len(sys.argv) < 5:
        print("usage: %s outfn num_bp disc_len L"%sys.argv[0])
        sys.exit()

    outfn        = sys.argv[1]

    nbp          = int(sys.argv[2])
    disc_len     = float(sys.argv[3])
    L            = float(sys.argv[4])

    scale_factor = 1./3.4

    ev_size = disc_len
    # ev_size = 0
    periodic_box = np.array([[0,L],[0,L],[0,L]])

    conf = gdna.gen_DNA_conf(nbp, disc_len, periodic_box=periodic_box, lb=40, ev_size=ev_size, first_pos=None, first_triad=None)

    scale_conf = conf*scale_factor
    scale_box  = periodic_box*scale_factor

    gen_conffile_with_DNA(
        outfn, 4, 1, 1, scale_box, scale_conf, [1], dna_bond_type = 1, dna_angle_type = 1)


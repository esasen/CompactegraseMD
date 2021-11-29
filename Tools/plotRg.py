#!/bin/env python3

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc

sys.path.append('inout')
sys.path.append('evals')
sys.path.append('conf')

import ReadXYZ as rxyz
import periodic_boundary as bp
    
    
    
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: %s fn" % sys.argv[0])
        sys.exit()

    fn   = sys.argv[1]

    # data = rxyz.LoadXYZ(fn,savenpy=False,loadnpy=False)
    data = rxyz.LoadXYZ(fn)
    dna_conf = rxyz.get_atoms_of_type(data, ['1','2'])

    print(np.shape(data['data']))
    print(np.shape(dna_conf))

    dna_uconf = bp.unwrap_dna(dna_conf)



    # uconf = bp.unwrap
    #
    #
    # if mmax < 1:
    #     mmax = 1
    #
    # conf = rxyz.ReadXYZ(fn)
    # lbs = clb.cal_persistence_length(conf,m_max=mmax)
    # print(lbs)


    

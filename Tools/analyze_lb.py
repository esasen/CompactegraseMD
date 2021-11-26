#!/bin/env python3

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import rc
   
import ReadXYZ as rxyz
import cal_persistence_length as clb
    
    
    
if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: %s fn mmax" % sys.argv[0])
        sys.exit()

    fn   = sys.argv[1]
    mmax = int(sys.argv[2])
    
    if mmax < 1:
        mmax = 1
    
    conf = rxyz.ReadXYZ(fn)
    lbs = clb.cal_persistence_length(conf,m_max=mmax)
    print(lbs)


    

import numpy as np
import sys,os
from FileRead import FileRead

"""
########################################################
Read/ReadXYZ.py
    
    specs,snapshots = LoadXYZ(filename)
    specs,snapshots = ReadXYZ(filename)
    
        ReadXYZ always reads the xyz file while LoadXYZ accesses the 
        binary if it exists and creates it if it doesn't such that access
        will be accelerated next time.
        
On Execution: 
    Read/ReadXYZ.py filename 
        
        Reads state file, prints specs in terminal and creates trajectory
        binary file.
        
########################################################
"""

def LoadXYZ(filename,savenpy=True,loadnpy=True):
    fnpy = '.'.join(filename.split('.')[:-1])+'.npy'
    if os.path.isfile(fnpy) and loadnpy:
        data = np.load(fnpy,allow_pickle=True)
    else:
        data = ReadXYZ(filename)
        if savenpy:
            SaveXYZ(fnpy,data)
    dat = dict()
    dat['data']  = data
    dat['types'] = ReadXYZ_atomtypes(filename)
    return dat

def ReadXYZ(fn):
    data = list()
    F = FileRead(fn)
    line = F.readline()
    while line!='':
        ll = F.linelist()
        if len(ll)>=4 and ll[0]!='Atoms.':
            snapshot = list()
            while len(ll)>=4:
                snapshot.append( [float(ft) for ft in ll[1:4]])
                line = F.readline()
                ll   = F.linelist()
            data.append(snapshot)
        line = F.readline()
    data = np.array(data)
    return data
    
def ReadXYZ_atomtypes(fn):
    data = list()
    F = FileRead(fn)
    line = F.readline()
    num = 0
    types = list()
    while line!='':
        ll = F.linelist()
        if len(ll)>=4 and ll[0]!='Atoms.':
            num += 1
            if num > 1:
                break
            while len(ll)>=4:
                types.append(ll[0])
                line = F.readline()
                ll   = F.linelist()
        line = F.readline()
    return types

def SaveXYZ(outname,data):
    if outname[-4:] == '.npy':
        outn = outname
    else:
        outn = outname + '.npy'
    np.save(outn,data)

def get_atoms_of_type(data,typeids):
    select = np.array([i for i in range(len(data['types'])) if data['types'][i] in typeids])
    return data['data'][:,select]

if __name__ == "__main__":
    
    if len(sys.argv) < 2:
        print("usage: python %s filename"%sys.argv[0])
        sys.exit(0)
    fn  = sys.argv[1]
    data = ReadXYZ(fn)
    fnpy = fn[:-4]+'.npy'
    SaveXYZ(fnpy,data)
    print(np.shape(data))

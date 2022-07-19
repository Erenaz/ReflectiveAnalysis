import numpy as np
import h5py
import sys, os
import json
import scipy.signal as sig
import ROOT
import root_numpy

import argparse
import matplotlib.pyplot as plt

import pickle as pkl

def getRefTemp(tempFn,nChan=4,nSamp=256,Eang=30.0,Hang=30.0,Cang=0.0):
    treenm="Templates"
    tree=ROOT.TChain(treenm)
    print("loading the Tree")
    tree.Add(tempFn)
    NEntries = tree.GetEntries()
    MetaArray = root_numpy.tree2array(tree,['EAng','HAng','coneAng'])

    DataArray = root_numpy.tree2array(tree,'wave.fData')
    DataArray = np.concatenate(DataArray)
    DataArray = np.reshape(DataArray,[NEntries,nChan,nSamp])
    DataArray = DataArray[:,1,:]

    iRef=None
    for i, angs in enumerate(MetaArray):
        #print('angs',float(angs[0]),float(angs[1]),float(angs[2]))
        #print('ref',float(Eang),float(Hang),float(Cang))
        if ((float(angs[0]) == float(Eang)) and (float(angs[1]) == float(Hang)) and (float(angs[2]) == float(Cang))):
            iRef = i

    if iRef:
        return DataArray[iRef]
    else:
        print('Error. Did not find reference template')
        sys.exit()






parser = argparse.ArgumentParser(description='Convert a root template into pickle format')
parser.add_argument('files', type=str, default='', help='File to run on')
parser.add_argument('saveAs', type=str)

args = parser.parse_args()
filesToRead = args.files
saveAs = args.saveAs



refTemp = getRefTemp(filesToRead)

plt.plot(refTemp)
plt.show()


pklTemplate = {}


#Dictionary format for templates is data[zenith][azimuth][viewing_angle]
#In this case, template provided is 30deg in E/Hplane
#So conversion is 30deg off of 90deg for both planes in radians
zen = np.deg2rad(120.0)
azi = np.deg2rad(30.0)
viewing_angle = np.deg2rad(0.0)
#pklTemplate[zen][azi][viewing_angle] = refTemp
pklTemplate[zen] = {}
pklTemplate[zen][azi] = {}
pklTemplate[zen][azi][viewing_angle] = refTemp

print(f'pkl is {pklTemplate}')

print(f'Saving file as {saveAs}')

with open(saveAs, 'wb') as fout:
    pkl.dump(pklTemplate, fout)

print(f'Saved file as {saveAs}')

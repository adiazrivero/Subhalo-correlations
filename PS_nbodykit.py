from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from nbodykit.source.catalog import HDFCatalog,CSVCatalog
from nbodykit import transform
from nbodykit.lab import ProjectedFFTPower, FFTPower

f = CSVCatalog('catalog_parttype1_*.txt', ['x','y','z'])

print f

f['Position'] =  transform.StackColumns(f['x'], f['y'], f['z'])
f['Weight'] = f['Weight'] * 2.8 * 1e4

f['r'] = np.linalg.norm(f['Position'],axis=1)
med = np.median(f['r'])
select = f['r'] == med

shift = f['Position'][select]
f['Position'] =  f['Position'] - shift

f['r'] = np.linalg.norm(f['Position'],axis=1)

bs = 600

select2 = f['r'] <= 600/2

subf = f[select2]

print subf

fig,ax = plt.subplots(1)
mesh = subf.to_mesh(Nmesh=128,BoxSize=bs,weight=f['Weight'])
k = ax.imshow(mesh.preview(axes=[0,1], Nmesh=128))
fig.colorbar(k)
#plt.savefig('mesh.png')
plt.show()


Pk = ProjectedFFTPower(mesh,BoxSize=bs, axes=(0,1), Nmesh=128)
Pk.save("projectfftpower.json")


"""
# plotting the power spectrum
r2 = ProjectedFFTPower.load("projectfftpower.json")
plt.loglog(r2.power['k'],r2.power['power'])
plt.show()
"""

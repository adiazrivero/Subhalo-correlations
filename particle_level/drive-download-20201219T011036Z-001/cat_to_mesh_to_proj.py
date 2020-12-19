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
from nbodykit.lab import ProjectedFFTPower, FFTPower,BigFileMesh
from nbodykit.source.mesh.array import ArrayMesh
from rotation import rotation
import dask.array as da
import cPickle as pickle

start_time = time.time()

inputs = sys.argv[1:]
for i in inputs:
    exec(i)

if numb == 78 or numb == 95:
    numb = '0' + str(numb)
else:
    numb = str(numb)

if name == 0:
    name = 'CDM'
elif name == 1:
    name = 'ETHOS_1'
elif name == 2:
    name = 'ETHOS_4'

print name

if numb == '127':
    print "z = 0"
elif numb == '095':
    print 'z = 0.5'
elif numb == '078':
    print 'z = 1'
else:
    print "wrong redshift specified"
    sys.exit()

lim = 300
bs = 2 * lim
nmesh = 1024
h = 0.6909
#mpart = 1.90382332e4 / h #1.90382332e-06 in the file['Header'].attrs['MassTable'][1]

dicti = pickle.load(open('/n/home04/adiazrivero/conv_powerspec/catalog_level/host_%s_%s.txt' % (name,numb),'rb'))
shiftx = dicti['host_center'][0]
shifty = dicti['host_center'][1]
shiftz = dicti['host_center'][2]

f = CSVCatalog('/n/dvorkin_lab/anadr/%s_catalog/%s/%s_%s_parttype1_*.txt' % (name,numb,name,numb), ['x','y','z'])
f['Position'] = transform.StackColumns(f['x'], f['y'],f['z'])
#f['Value'] = mpart

N_part = len(f)
print 'number of particles: %s' % N_part

box_vol = bs**3
nbar = N_part / box_vol

print nbar

f['x'] = f['x'] / h - shiftx
f['y'] = f['y'] / h - shifty
f['z'] = f['z'] / h - shiftz

select1 = f['x'] < lim
select2 = f['x'] > -lim
select3 = f['y'] < lim
select4 = f['y'] > -lim
select5 = f['z'] < lim
select6 = f['z'] > -lim

subcat = f[select1 & select2 & select3 & select4 & select5 & select6]
subcat['Position'] =  transform.StackColumns(subcat['x'], subcat['y'], subcat['z'])

"""print subcat 

N_part = len(subcat)
print 'number of particles: %s' % N_part

box_vol = bs**3
nbar = N_part / box_vol

print nbar

dicti = {'nbar': nbar, 'N_part': N_part}

with open('/n/home04/adiazrivero/conv_powerspec/particle_level/nbar/nbar_%s_%s.txt' % (name,numb), 'w') as file:
     file.write(pickle.dumps(dicti))

#mesh = subcat.to_mesh(Nmesh=nmesh,BoxSize=bs)
#mesh.save('/n/dvorkin_lab/anadr/%s_%s_%s_%s_hbugfix.bigfile' % (name,numb,nmesh,bs,h))
#mesh.save('/n/dvorkin_lab/anadr/%s_%s_%s_%s_h%s_mass.bigfile' % (name,numb,nmesh,bs,h))
#mesh.save('/n/dvorkin_lab/anadr/%s_%s_%s_%s_h%s.bigfile' % (name,numb,nmesh,bs,h))"""

mesh = BigFileMesh('/n/dvorkin_lab/anadr/%s_%s_%s_%s_hbugfix.bigfile' % (name,numb,nmesh,bs),'Field')

print("%s seconds" % (time.time() - start_time))

count = 0

while count < num_proj:
   count += 1
   nnx = np.random.uniform(0,10)
   nny = np.random.uniform(0,10)
   nnz = np.random.uniform(0,10)
   theta = np.random.uniform(0,2*np.pi)

   nx = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnx
   ny = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nny
   nz = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnz

   R = rotation(nx,ny,nz,theta)

   subcat['Position'] = da.transpose(da.dot(R,da.transpose(subcat['Position'])))

   mesh = subcat.to_mesh(Nmesh=nmesh,BoxSize=bs)

   proj1 = np.fft.fftshift(mesh.preview(axes=[0,1], Nmesh=nmesh))
   proj2 = np.fft.fftshift(mesh.preview(axes=[0,2], Nmesh=nmesh))
   proj3 = np.fft.fftshift(mesh.preview(axes=[1,2], Nmesh=nmesh))

   np.save('/n/dvorkin_lab/anadr/%s_proj/%s/hbugfix/%s_%s_projxy_%s_%s' % (name,numb,name,numb,count,item),proj1)
   np.save('/n/dvorkin_lab/anadr/%s_proj/%s/hbugfix/%s_%s_projxz_%s_%s' % (name,numb,name,numb,count,item),proj2)
   np.save('/n/dvorkin_lab/anadr/%s_proj/%s/hbugfix/%s_%s_projyz_%s_%s' % (name,numb,name,numb,count,item),proj3)

   print("%s seconds" % (time.time() - start_time))

print("%s seconds" % (time.time() - start_time))

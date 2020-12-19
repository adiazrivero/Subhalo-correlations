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
import dask.array as da
from rotation import *
from scipy.ndimage import rotate

#plt.rc('text', usetex=True)

name = 'CDM'
numb = '095'
bs = 600
nmesh = 1024

num = 200

#mesh = BigFileMesh('ethos4_mesh_%s_%s_h07.bigfile' % (nmesh,bs),'Field')
mesh1 = BigFileMesh('/n/dvorkin_lab/anadr/%s_%s_1024_600_h0.6909.bigfile' % (name,numb),'Field')
#mesh2 = BigFileMesh('/n/dvorkin_lab/anadr/%s_%s_1024_600_hbugfix.bigfile' % (name,numb),'Field')

mesh1 = np.fft.fftshift(mesh1.preview(Nmesh=nmesh))
print 'done'

rot = rotate(mesh1, 30, reshape=False)

sys.exit()

nnx = np.random.uniform(0,10)
nny = np.random.uniform(0,10)
nnz = np.random.uniform(0,10)
theta = np.random.uniform(0,2*np.pi)

nx = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnx
ny = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nny
nz = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnz

R = rotation(nx,ny,nz,theta)

mesh2 = np.dot(R,mesh1)

sys.exit()

#mesh2 = da.transpose(da.dot(R,da.transpose(subcat['Position'])))
mesh2 = da.transpose(da.dot(R,da.transpose(mesh1)))


sys.exit()

proj1 = np.fft.fftshift(mesh1.preview(axes=[0,1], Nmesh=nmesh))
proj1 = proj1[num:-num,num:-num]


# Generate a random projection angle
theta = np.random.uniform(0, np.pi, size=num_maps)
phi = np.random.uniform(0, 2*np.pi, size=num_maps)

theta_hat = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi),-np.sin(theta)]).T
phi_hat = np.array([-np.sin(phi),np.cos(phi),np.zeros(num_maps)]).T
    
# Compute projected positions of subhalos
sub_pos = np.vstack((np.tensordot(pos,theta_hat[imap],axes=1),np.tensordot(pos,phi_hat[imap],axes=1))).T

#proj2 = np.fft.fftshift(mesh2.preview(axes=[0,1], Nmesh=nmesh))

fig,ax = plt.subplots(1)
k = ax.imshow(np.log10(proj1),extent=[-bs/2,bs/2,-bs/2,bs/2],cmap='inferno')
fig.colorbar(k)
#plt.show()
plt.savefig('mesh_old.png')

sys.exit()

fig,ax = plt.subplots(1)
k = ax.imshow(np.log10(proj2),extent=[-bs/2,bs/2,-bs/2,bs/2],cmap='inferno')
fig.colorbar(k)
#plt.show()
plt.savefig('mesh_new.png')

sys.exit()

name = 'ETHOS_4'
numb = '095'

mesh = BigFileMesh('/n/dvorkin_lab/anadr/%s_%s_1024_600_h0.6909.bigfile' % (name,numb),'Field')

proj = np.fft.fftshift(mesh.preview(axes=[0,1], Nmesh=nmesh))

select = proj <= 0
proj[select] = 1e-4

fig,ax = plt.subplots(1)
k = ax.imshow(np.log10(proj),extent=[-bs/2,bs/2,-bs/2,bs/2],cmap='inferno')
#fig.colorbar(k)
#plt.show()
plt.savefig('ETHOS_4_part.png')

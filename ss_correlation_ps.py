from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import h5py
import sys
import time
from test import test

start_time = time.time()

list = ['fof_subhalo_tab_127.0.hdf5','fof_subhalo_tab_127.1.hdf5','fof_subhalo_tab_127.2.hdf5','fof_subhalo_tab_127.3.hdf5','fof_subhalo_tab_127.4.hdf5','fof_subhalo_tab_127.5.hdf5','fof_subhalo_tab_127.6.hdf5','fof_subhalo_tab_127.7.hdf5','fof_subhalo_tab_127.8.hdf5','fof_subhalo_tab_127.9.hdf5','fof_subhalo_tab_127.10.hdf5','fof_subhalo_tab_127.11.hdf5','fof_subhalo_tab_127.12.hdf5','fof_subhalo_tab_127.13.hdf5','fof_subhalo_tab_127.14.hdf5','fof_subhalo_tab_127.15.hdf5']

masses = []
positions1 = []
vmax1 = []
for i in list:
    file = h5py.File(i,'r')
    if len(file['Subhalo'].keys()) != 0:
        masses.append(file['Subhalo']['SubhaloMass'].value) #in 10^10 M_sun/h
        positions1.append(file['Subhalo']['SubhaloPos'].value) # in kpc/h, box coordinates

masses = [i for j in masses for i in j]
positions1 = [i for j in positions1 for i in j]
parent_mass = masses[0]
parent_pos = positions1[0]
positions = [i - parent_pos for i in positions1]

subh = zip(masses[1:],positions[1:])

Rnge = 600
subh_cut=[]
[subh_cut.append(i) for i in subh if i[0] > 1.5e-4 and -Rnge/2. < i[1][0] < Rnge/2. and -Rnge/2. < i[1][1] < Rnge/2. and -Rnge/2. < i[1][2] < Rnge/2.]

subh_mass,subh_pos = zip(*subh_cut)
#print "total number of ETHOS subhalos within r < %s kpc/h: %s " % (Rnge/2,len(subh_cut))

def rotation(nx,ny,nz,theta):
    R = [[np.cos(theta) + (nx**2)*(1-np.cos(theta)) , nx*ny*(1-np.cos(theta)) - nz*np.sin(theta) , nx*nz*(1-np.cos(theta)) + ny*np.sin(theta)], [nx*ny*(1-np.cos(theta)) + nz*np.sin(theta) , np.cos(theta) + (ny**2)*(1-np.cos(theta)) , ny*nz*(1-np.cos(theta)) - nx*np.sin(theta)], [nz*nx*(1-np.cos(theta)) - ny*np.sin(theta) , nz*ny*(1-np.cos(theta)) + nx*np.sin(theta) , np.cos(theta) + (nz**2)*(1-np.cos(theta))]]

    return R

rnge = 100
coords = []
count = 0
while count < 1000:
    count += 1
    nnx = np.random.uniform(0,10)
    nny = np.random.uniform(0,10)
    nnz = np.random.uniform(0,10)
    theta = np.random.uniform(0,2*np.pi)

    nx = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnx
    ny = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nny
    nz = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnz
    #print "vector magnitude: %s" % (nx**2+ny**2+nz**2) #verifying it is a unit vector

    R = rotation(nx,ny,nz,theta)
    #print np.linalg.det(R) #verifying det=1
    rot_pos = [np.dot(R,i) for i in subh_pos]

    proj_xy = [[i[0],i[1]] for i in rot_pos if -rnge/2. < i[0] < rnge/2. and -rnge/2. < i[1] < rnge/2.]
    proj_xz = [[i[0],i[2]] for i in rot_pos if -rnge/2. < i[0] < rnge/2. and -rnge/2. < i[2] < rnge/2.]
    proj_yz = [[i[1],i[2]] for i in rot_pos if -rnge/2. < i[1] < rnge/2. and -rnge/2. < i[2] < rnge/2.]

    coords.append(proj_xy)
    coords.append(proj_xz)
    coords.append(proj_yz)

tot_num_subh = []
for i in coords:
    tot_num_subh.append(len(i))

#print "average number of subhalos within r < 50 kpc/h after rotating & projecting: %s" % np.mean(tot_num_subh)

bns = 20
bin_size = rnge/bns
bin_avg = []
coadd_ps = []
coadd_corr = []

for i in coords:
    x,y = zip(*i)
    H,xedges,yedges = np.histogram2d(x,y,bins=bns)
    k_val = [(2*np.pi)/i for i in xedges]
    #print k_val
    Nbar = H.mean()
    bin_avg.append(Nbar)
    corr = (H - Nbar) / Nbar
    coadd_corr.append(corr)
    ft = np.fft.fft2(corr)
    ft = [i/(2*np.pi) for i in ft] #numpy fft convention: e^(-i2\pixk)
    ps2D = np.abs(ft)**2
    coadd_ps.append(ps2D)

#print "average nbar: %s " % np.mean(bin_avg)

show_nr = False
show_ps = True

tot_corr = sum(coadd_corr)
tot_corr = [i/len(coadd_corr) for i in tot_corr]

if show_nr == True:
    py.figure('(n - nbar)/nbar')
    plt.imshow(tot_corr,extent=[-rnge/2, rnge/2, -rnge/2, rnge/2],interpolation='nearest')
    plt.colorbar()
    py.show()

coadd_ps = [i/len(coadd_ps) for i in coadd_ps]
tot_ps = sum(coadd_ps)
tot_ps = np.fft.fftshift(tot_ps)
kx = np.fft.fftfreq(tot_ps.shape[0],d=bin_size)
kx = [2*np.pi*i for i in kx]
ky = np.fft.fftfreq(tot_ps.shape[1],d=bin_size)
ky = [2*np.pi*i for i in ky]

if show_ps == True:
    py.figure('2d Power Spectrum')
    py.imshow(tot_ps,extent=[min(kx),max(kx),min(ky),max(ky)],interpolation='nearest')
    plt.colorbar()
    py.show()

def oneD_ps(data,pixel_size=1):
    ps2d = np.array(data)
    pixx, pixy = ps2d.shape
    x1 = np.arange(-pixx/2.,pixx/2.)
    y1 = np.arange(-pixy/2.,pixy/2.)
    x,y = np.meshgrid(y1,x1)

    k = np.sqrt(x**2+y**2)
    kmax = k.max()
    dk = pixel_size

    K1 = np.arange(kmax/dk)*dk + dk/2.
    conv = kx[1]/K1[-1]
    K = K1*conv
    nk = len(K)
    print K

    ring_mean = []
    for i in range(nk):
        kmin = i*dk
        kmax = kmin + dk
        ind = (k>=kmin) * (k<kmax)
        #print ind*1
        ring_mean.append(data[ind].mean())

    return ring_mean,K

pix_size_k = 2*np.pi/bin_size
ps1d,K = oneD_ps(tot_ps,pixel_size=pix_size_k)

A_pix = bin_size**2 #(kpc/h)^2
A_box = rnge**2 #(kpc/h)^2
norm = (A_pix**2)/(A_box)
ps1d = [norm*i for i in ps1d]
avg_ps = np.mean(ps1d)
var = [(i**2-avg_ps)**2/len(coords) for i in ps1d]

#tst = test()
tst = None

py.figure('1d Power Spectrum')
ax = plt.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
plt.errorbar(K, ps1d, xerr=0, yerr=var)
if tst is not None:
    plt.loglog(K,tst)
ax.set_xlim(min(K)-1e-5,max(K))
ax.set_xlabel('k [h kpc]^-1')
ax.set_ylabel('P_ss(k)[kpc/h]^2')
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

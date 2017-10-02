import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
import h5py
import sys
import time

#######################################################################
#importing data
#######################################################################

list = ['fof_subhalo_tab_127.0.hdf5','fof_subhalo_tab_127.1.hdf5','fof_subhalo_tab_127.2.hdf5','fof_subhalo_tab_127.3.hdf5','fof_subhalo_tab_127.4.hdf5','fof_subhalo_tab_127.5.hdf5','fof_subhalo_tab_127.6.hdf5','fof_subhalo_tab_127.7.hdf5','fof_subhalo_tab_127.8.hdf5','fof_subhalo_tab_127.9.hdf5','fof_subhalo_tab_127.10.hdf5','fof_subhalo_tab_127.11.hdf5','fof_subhalo_tab_127.12.hdf5','fof_subhalo_tab_127.13.hdf5','fof_subhalo_tab_127.14.hdf5','fof_subhalo_tab_127.15.hdf5']

masses = []
positions1 = []
vmax1 = []
for i in list:
    file = h5py.File(i,'r')
    #print "%s keys: %s" % (i,file['Subhalo'].keys())
    if len(file['Subhalo'].keys()) != 0:
        masses.append(file['Subhalo']['SubhaloMass'].value) #in 10^10 M_sun/h
        positions1.append(file['Subhalo']['SubhaloPos'].value) # in kpc/h, box coordinates
        vmax1.append(file['Subhalo']['SubhaloVmax'].value)

masses = [i for j in masses for i in j]
positions1 = [i for j in positions1 for i in j]
vmax = [i for j in vmax1 for i in j]
parent_mass = masses[0]
parent_pos = positions1[0]
positions = [i - parent_pos for i in positions1] #halocentric coordinates

subh = zip(masses[1:],positions[1:],vmax[1:])

subh_cut=[]
[subh_cut.append(i) for i in subh if i[0] > 1.5e-4 and i[1][0] < 300 and i[1][1] < 300 and i[1][2] < 300] #and i[0] < 1e-2]

subh_mass,subh_pos,subh_vmax = zip(*subh_cut)
#subh_mass = [i*1e10 for i in subh_mass]

"""counts,bins,bars = pl.hist(subh_mass,log=True,bins=np.logspace(4.0,11.0,150))
f = pl.gca()
f.set_title("Number of subhalos per mass bin in ETHOS")
f.set_xscale("log")
f.set_xlabel("Mass [$M_{\odot}$]")
f.set_ylabel("N")
pl.show()

counts2,bins2,bars2 = pl.hist(subh_vmax,log=True,bins=np.logspace(0,3.0,150))
f = pl.gca()
f.set_title("Number of subhalos per V_max bin in ETHOS")
f.set_xscale("log")
f.set_xlabel("Vmax")
f.set_ylabel("N")
pl.show()

#m = np.arange(10**(5),10**(11),10**5)
#dndm = (1e-5)*(m/(2.52*10**7))**(-1.9)

plt.loglog(bins[:149],counts)
#plt.plot(m,dndm,color='g',label='$dNdm \propto m^{-1.9}$')
plt.xlim(1e5,1e11)
plt.ylabel('dN/dM')
plt.xlabel('Mass [$M_{\odot}$]')
plt.show()

plt.loglog(bins2[:149],counts2)
#plt.plot(m,dndm,color='g',label='$dNdm \propto m^{-1.9}$')
plt.xlim(1,1e3)
plt.ylabel('dN/dVmax')
plt.xlabel('Vmax')
plt.show()"""

#######################################################################
#rotating,projecting,fft
#######################################################################

def rotation(nx,ny,nz,theta):
    R = [[np.cos(theta) + (nx**2)*(1-np.cos(theta)) , nx*ny*(1-np.cos(theta)) - nz*np.sin(theta) , nx*nz*(1-np.cos(theta)) + ny*np.sin(theta)], [nx*ny*(1-np.cos(theta)) + nz*np.sin(theta) , np.cos(theta) + (ny**2)*(1-np.cos(theta)) , ny*nz*(1-np.cos(theta)) - nx*np.sin(theta)], [nz*nx*(1-np.cos(theta)) - ny*np.sin(theta) , nz*ny*(1-np.cos(theta)) + nx*np.sin(theta) , np.cos(theta) + (nz**2)*(1-np.cos(theta))]]
    return R

start_time = time.time()

coords_10 = []
coords_20 = []
coords_30 = []
coords_50 = []
count = 0

while count < 1000:
    count += 1
    nnx = np.random.uniform(0,10)
    nny = np.random.uniform(0,10)
    nnz = np.random.uniform(0,10)
    theta = np.random.uniform(0,2*np.pi)

    #normalizing \vec{n}
    nx = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnx
    ny = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nny
    nz = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnz

    '''print 'unit vector: (%s,%s,%s) ' % (nx,ny,nz)
    print 'magnitude: %s' % np.sqrt(nx**2 + ny**2 + nz**2)'''

    R = rotation(nx,ny,nz,theta)
    rot_pos = [np.dot(R,i) for i in subh_pos]

    proj_xy = [[i[0],i[1]] for i in rot_pos if -10 < i[0] < 10 and -10 < i[1] < 10]
    proj_xz = [[i[0],i[2]] for i in rot_pos if -10 < i[0] < 10 and -10 < i[2] < 10]
    proj_yz = [[i[1],i[2]] for i in rot_pos if -10 < i[1] < 10 and -10 < i[2] < 10]

    coords_10.append(proj_xy)
    coords_10.append(proj_xz)
    coords_10.append(proj_yz)

    proj_xy_2 = [[i[0],i[1]] for i in rot_pos if -20 < i[0] < 20 and -20 < i[1] < 20]
    proj_xz_2 = [[i[0],i[2]] for i in rot_pos if -20 < i[0] < 20 and -20 < i[2] < 20]
    proj_yz_2 = [[i[1],i[2]] for i in rot_pos if -20 < i[1] < 20 and -20 < i[2] < 20]

    coords_20.append(proj_xy_2)
    coords_20.append(proj_xz_2)
    coords_20.append(proj_yz_2)

    proj_xy_3 = [[i[0],i[1]] for i in rot_pos if -30 < i[0] < 30 and -30 < i[1] < 30]
    proj_xz_3 = [[i[0],i[2]] for i in rot_pos if -30 < i[0] < 30 and -30 < i[2] < 30]
    proj_yz_3 = [[i[1],i[2]] for i in rot_pos if -30 < i[1] < 30 and -30 < i[2] < 30]

    coords_30.append(proj_xy_3)
    coords_30.append(proj_xz_3)
    coords_30.append(proj_yz_3)

    proj_xy_5 = [[i[0],i[1]] for i in rot_pos if -50 < i[0] < 50 and -50 < i[1] < 50]
    proj_xz_5 = [[i[0],i[2]] for i in rot_pos if -50 < i[0] < 50 and -50 < i[2] < 50]
    proj_yz_5 = [[i[1],i[2]] for i in rot_pos if -50 < i[1] < 50 and -50 < i[2] < 50]

    coords_50.append(proj_xy_5)
    coords_50.append(proj_xz_5)
    coords_50.append(proj_yz_5)

print("--- %s seconds ---" % (time.time() - start_time))

file = open('coords_10x10_2.py', 'w')
file.write('%s\n' % coords_10)
file.close()

file = open('coords_20x20_2.py', 'w')
file.write('%s\n' % coords_20)
file.close()

file = open('coords_30x30_2.py', 'w')
file.write('%s\n' % coords_30)
file.close()

file = open('coords_50x50_2.py', 'w')
file.write('%s\n' % coords_50)
file.close()

print("--- %s seconds ---" % (time.time() - start_time))

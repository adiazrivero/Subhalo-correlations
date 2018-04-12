from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from R_ein import rein_sigmac
from functions import simul_data,rotation,projections,angular_average,twoD_ps

start_time = time.time()

# importing data from halo catalogs

list1 = []
for i in range(16):
    #list1.append('fof_subhalo_tab_127.%s.hdf5' % i)
    list1.append('/Users/anadiazrivero/Desktop/ETHOS/CDM/fof_subhalo_tab_127.%s.hdf5' % i)

shft = 0
rnge = 100
num_proj = 10
cut = True

mass,pos,rh,rvmax = simul_data(list1,mhigh_cut=cut,Rnge=600)

print "total number of subhalos within r_3d < 300 kpc: %s " % np.shape(mass)

mass2 = [i*10**(10) for i in mass] #units of M_sun

print "max mass = %s" % max(mass2)
print "min mass = %s" %  min(mass2)

Mhost = 1.04e12
zs = 1
zl = 0.5
rein,sigmac = rein_sigmac(Mhost,zs,zl) #in kpc and M_sun/kpc^2

print "r_ein = %s " % rein
print "sigma_crit = %s " % sigmac

posmass,_ = projections(pos,mass2,rh,rvmax,rnge=rnge,shift=shft,num_proj=num_proj)

posmass2 = np.asarray(posmass)
print "proj num = %s " % len(posmass2)

np.save('projections',posmass2)

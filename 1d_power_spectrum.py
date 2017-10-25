from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
#from test import test
from functions import simul_data,rotation,projections,twoD_ps,oneD_ps,error_bars

start_time = time.time()

#tst = test()
tst = None

subh_mass,subh_pos = simul_data(Rnge=600)

coords = projections(subh_pos,rnge=100,shift=0,num_proj=1000)

#coadd_ps,tot_ps,kx,ky,bin_size = twoD_ps(coords,bns=20,rnge=100,shift=0,show_nr=False,show_ps=False)
coadd_ps,tot_ps,kx,ky,bin_size = twoD_ps(coords,bns=20,rnge=100,shift=0,show_nr=False,show_ps=False)


pix_size_k = np.abs(kx[0]-kx[1])
ps1d,K,norm = oneD_ps(tot_ps,kx,ky,rnge=100,bin_size=bin_size,pixel_size=pix_size_k)

error = error_bars(coadd_ps,ps1d,kx,ky,rnge=100,pix_size=pix_size_k)

py.figure('1d Power Spectrum')
ax = plt.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
"""if error == False:
    plt.loglog(K,ps1d)"""
plt.errorbar(K, ps1d, yerr=error)
if tst is not None:
    plt.loglog(K,tst)
#ax.set_xlim(1e-3,max(K)+1e-1)
ax.set_ylim(1e-1,max(ps1d)+1)
ax.set_xlabel('k [h kpc]^-1')
ax.set_ylabel('P_ss(k)[kpc/h]^2')
plt.show()


print("--- %s seconds ---" % (time.time() - start_time))

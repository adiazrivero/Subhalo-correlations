from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
from functions import simul_data,rotation,projections,twoD_ps,angular_average,variance,poisson_realization

start_time = time.time()

rnge = 100
list1 = ['fof_subhalo_tab_127.0.hdf5','fof_subhalo_tab_127.1.hdf5','fof_subhalo_tab_127.2.hdf5','fof_subhalo_tab_127.3.hdf5','fof_subhalo_tab_127.4.hdf5','fof_subhalo_tab_127.5.hdf5','fof_subhalo_tab_127.6.hdf5','fof_subhalo_tab_127.7.hdf5','fof_subhalo_tab_127.8.hdf5','fof_subhalo_tab_127.9.hdf5','fof_subhalo_tab_127.10.hdf5','fof_subhalo_tab_127.11.hdf5','fof_subhalo_tab_127.12.hdf5','fof_subhalo_tab_127.13.hdf5','fof_subhalo_tab_127.14.hdf5','fof_subhalo_tab_127.15.hdf5']

# SIMULATION DATA

data = True
if data == True:
    subh_mass,subh_pos = simul_data(list1,Rnge=600)
    coords,avg_num_subh = projections(subh_pos,rnge=rnge,shift=0,num_proj=1000)
    individual_ps,tot_ps,kx,ky,bin_size = twoD_ps(coords,bns=21,rnge=rnge,shift=0,show_nr=False,show_ps=False)
    pix_size_k = np.abs(kx[0]-kx[1])
    ps1d,K,norm = angular_average(tot_ps,kx,ky,rnge=rnge,bin_size=bin_size,dk=pix_size_k)
    var = variance(individual_ps,ps1d,len(coords),kx,ky,rnge=rnge,pix_size=pix_size_k)
    std = [np.sqrt(i) for i in var]
    #std = None

    print "P_ss(k) = %s " % ps1d
    print "var= %s" % var
    print "std= %s" % std

# POISSON-DISTRIBUTED REALIZATIONS

test = True
if test == True:
    """vals = [0.3,0.5,1,2,3]
    for x in vals:
        print "_____________________________"
        print x
        num_subh = avg_num_subh*x
        print "average number of subhalos on the plane (poisson): %s" % num_subh
        poiss_coords = poisson_realization(rnge=rnge,subh_num=num_subh,num_proj=3000)
        individual_ps2,tot_ps2,kx2,ky2,bin_size2 = twoD_ps(poiss_coords,bns=21,rnge=rnge,shift=0,show_nr=False,show_ps=False)

        pix_size_k2 = np.abs(kx[0]-kx[1])
        noise,K2,norm2 = angular_average(tot_ps2,kx2,ky2,rnge=rnge,bin_size=bin_size,dk=pix_size_k2)

        print noise
        var2 = variance(individual_ps2,noise,len(poiss_coords),kx2,ky2,rnge=rnge,pix_size=pix_size_k2)
        std2 = [np.sqrt(i) for i in var]
        #std2 = None
        print "var= %s" % var2
        print "std= %s" % std2"""

    num_subh = avg_num_subh
    print "average number of subhalos on the plane (poisson): %s" % num_subh
    poiss_coords = poisson_realization(rnge=rnge,subh_num=num_subh,num_proj=3000)
    individual_ps2,tot_ps2,kx2,ky2,bin_size2 = twoD_ps(poiss_coords,bns=21,rnge=rnge,shift=0,show_nr=False,show_ps=False)
    pix_size_k2 = np.abs(kx[0]-kx[1])
    noise,K2,norm2 = angular_average(tot_ps2,kx2,ky2,rnge=rnge,bin_size=bin_size,dk=pix_size_k2)
    var2 = variance(individual_ps2,noise,len(poiss_coords),kx2,ky2,rnge=rnge,pix_size=pix_size_k2)
    std2 = [np.sqrt(i) for i in var2]
    #std2 = None

    print "noise PS = %s " % noise
    print "var = %s" % var2
    print "std = %s" % std2

py.figure('1d Power Spectrum')
ax = plt.subplot(111)
ax.set_xscale("log")
ax.set_yscale("log")
if data == True:
    if std is None:
        plt.plot(K,ps1d,label='P_ss(k)')
    else:
        plt.errorbar(K, ps1d, yerr=std,label='P_ss(k)')
if test == True:
    if std2 is None:
        plt.plot(K2,noise,label='P_noise(k)')
    else:
        plt.errorbar(K2, noise, yerr=std2,label='P_noise(k)')
ax.set_xlabel('k [h kpc]^-1')
ax.set_ylabel('P(k)[kpc/h]^2')
plt.legend()
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

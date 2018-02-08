from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
from functions import simul_data,rotation,projections,twoD_nr,twoD_ps,angular_average,variance,poisson_realization,interpolation

np.set_printoptions(suppress=True,linewidth=np.nan,threshold=np.nan)

start_time = time.time()

shft = 0
rnge = 100
pix_num = 101
pix_size = rnge/pix_num

list1 = ['fof_subhalo_tab_127.0.hdf5','fof_subhalo_tab_127.1.hdf5','fof_subhalo_tab_127.2.hdf5','fof_subhalo_tab_127.3.hdf5','fof_subhalo_tab_127.4.hdf5','fof_subhalo_tab_127.5.hdf5','fof_subhalo_tab_127.6.hdf5','fof_subhalo_tab_127.7.hdf5','fof_subhalo_tab_127.8.hdf5','fof_subhalo_tab_127.9.hdf5','fof_subhalo_tab_127.10.hdf5','fof_subhalo_tab_127.11.hdf5','fof_subhalo_tab_127.12.hdf5','fof_subhalo_tab_127.13.hdf5','fof_subhalo_tab_127.14.hdf5','fof_subhalo_tab_127.15.hdf5']

data = True
test = True

std = None
std2 = None
plot_pss = True

# SIMULATION DATA
if data == True:
    subh_mass,subh_pos,rh,rvmax = simul_data(list1,Rnge=600,mhigh_cut=False)
    coordsm,avg_num_subh = projections(subh_pos,subh_mass,rnge=rnge,shift=shft,num_proj=1000)

    coords = []
    for i in coordsm:
        r,m = zip(*i)
        coords.append(r)

    ind_corr,coadd_corr,nbar = twoD_nr(coords,pix_num=pix_num,rnge=rnge,shift=shft,show_nr=False)
    individual_ps,tot_ps,kx,ky = twoD_ps(data=None,ind_data=ind_corr,pix_size=pix_size,rnge=rnge,shift=shft,show_ps=False)
    pix_size_k = np.abs(kx[0]-kx[1])
    ps1d,K = angular_average(tot_ps,kx,ky,rnge=rnge,pix_num=pix_num,dr=pix_size_k,remove_first=True)
    var = variance(individual_ps,ps1d,len(coords),kx,ky,rnge=rnge,pix_num=pix_num,pix_size=pix_size_k)
    std = [np.sqrt(i) for i in var]

# RANDOM REALIZATIONS
if test == True:
    if data == False:
        num_subh = 1406
    if data == True:
        num_subh = avg_num_subh

    poiss_coords = poisson_realization(rnge=rnge,subh_num=num_subh,num_proj=3000)
    ind_corr_poiss,tot_corr_poiss,nbar2 = twoD_nr(poiss_coords,pix_num=pix_num,rnge=rnge,shift=shft,show_nr=False)

    ind_ps2,tot_ps2,kx2,ky2 = twoD_ps(data=None,ind_data=ind_corr_poiss,pix_size=pix_size,rnge=rnge,shift=shft,show_ps=False)
    pix_size_k2 = np.abs(kx2[0]-kx2[1])
    noise,K2 = angular_average(tot_ps2,kx2,ky2,rnge=rnge,pix_num=pix_num,dr=pix_size_k2,remove_first=True)
    var2 = variance(ind_ps2,noise,len(poiss_coords),kx2,ky2,rnge=rnge,pix_num=pix_num,pix_size=pix_size_k2)
    std2 = [np.sqrt(i) for i in var2]

if plot_pss == True:
    py.figure('1d Power Spectrum (shift %s kpc/h), %s x %s pixels' % (shft,pix_num,pix_num))
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
    ax.set_xlim(4e-2,1)
    ax.set_xlabel('k [h kpc]^-1')
    ax.set_ylabel('P(k)[kpc/h]^2')
    if test == True:
        plt.axhline(y=1/nbar2,color='r',linestyle='-',label='1/nbar')
    plt.legend()
    plt.show()
    plt.clf()

export_p2sh = True

if export_p2sh == True:
    p2sh = [i*0.0048**2 for i in ps1d]
    file_p2sh = open('p2sh_sim.txt', 'w')
    for i in p2sh:
        file_p2sh.write('%s\n' % i)
    file_p2sh.close()

    file_k = open('k_sim.txt', 'w')
    for i in K:
        file_k.write('%s\n' % i)
    file_k.close()

print("--- %s seconds ---" % (time.time() - start_time))

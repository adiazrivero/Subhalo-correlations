from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pylab as py
import sys
import time
import h5py
import cmath
from rein import rein_sigmac
from functions import simul_data,rotation,projections,twoD_ps,angular_average

start_time = time.time()

rnge = 100
shft = 0


# importing data from halo catalogs
list1 = ['fof_subhalo_tab_127.0.hdf5','fof_subhalo_tab_127.1.hdf5','fof_subhalo_tab_127.2.hdf5','fof_subhalo_tab_127.3.hdf5','fof_subhalo_tab_127.4.hdf5','fof_subhalo_tab_127.5.hdf5','fof_subhalo_tab_127.6.hdf5','fof_subhalo_tab_127.7.hdf5','fof_subhalo_tab_127.8.hdf5','fof_subhalo_tab_127.9.hdf5','fof_subhalo_tab_127.10.hdf5','fof_subhalo_tab_127.11.hdf5','fof_subhalo_tab_127.12.hdf5','fof_subhalo_tab_127.13.hdf5','fof_subhalo_tab_127.14.hdf5','fof_subhalo_tab_127.15.hdf5']

mass,pos,rh,rvmax = simul_data(list1,Rnge=600)
mass2 = [i*10**(10) for i in mass] #units of M_sun
rein,sigmac = rein_sigmac(1.04e12,1,0.5) #in kpc and M_sun/kpc^2
h = 0.7
rein_kpch = h*rein
num_proj = 30

posmass,num_sub = projections(pos,mass2,rnge=rnge,shift=shft,num_proj=num_proj)

pos_2d = []
for i in posmass:
    r,m = zip(*i)
    pos_2d.append(r)

#########################################################################
#kappa avg using eq. 14 in our paper \kappa = (Nsub*mavg)/(A*sigma_c)
#########################################################################

arr = []
mav = []
for b in posmass:
    lst = [a for a in b if np.sqrt(a[0][0]**2+a[0][1]**2) <= rein_kpch]
    arr.append(lst)
    m = [a[1] for a in lst]
    mav.append(m)

mavg = np.asarray([np.mean(i) for i in mav])
N = np.asarray([len(i) for i in arr])

def kappa_avg(nsub,mavg,area,sigmac):
    kavg = (nsub*mavg)/(area*sigmac)
    return kavg

A = np.pi*rein_kpch**2
kavg = kappa_avg(N,mavg,A,sigmac)
print "average kappa in SL region (14): %s " % np.mean(kavg)

#####################################
#####################################

def tau(rs,rt):
    """
    takes in arrays of scale and tidal radii and returns array of tau
    """
    tau = [j/i for i,j in zip(rs,rt)]
    return tau

def mnfw(m,tau):
    """
    taken from Baltz et al (2009)
    """
    n = [(i**2/(i**2+1)**2) * ((i**2-1)*np.log(i)+i*np.pi-(i**2+1)) for i in tau]
    mnfw = [i/j for i,j in zip(m,n)]
    return mnfw

def conv_tnfw(x,y,tau,mnfw,rs):
    """
    inputs:
        x,y: x and y positions at which the convergence is evaluated
        tau: see above
        mnfw: see above
        rs: scale radius
    output:
        convergence for tNFW subhalo
    """
    z = np.sqrt(x**2+y**2)/rs
    Fz = []
    Lz = []
    for i in z:
        a = [-cmath.acos(1/j)/cmath.sqrt(j**2-1) if j < 1 else cmath.acos(1/j)/cmath.sqrt(j**2-1) for j in i]
        b = [cmath.log(j/(tau+np.sqrt(tau**2+j**2))) for j in i]
        Fz.append(a)
        Lz.append(b)
    Fz = np.asarray(Fz)
    Lz = np.asarray(Lz)
    sig = (mnfw/(sigmac*rs**2))*(tau**2/(2*np.pi*(tau**2+1)**2))*(((tau**2+1)/(z**2-1))*(1-Fz)+2*Fz-np.pi/np.sqrt(tau**2+z**2)+((tau**2-1)/(tau*np.sqrt(tau**2+z**2)))*Lz)
    return sig

v = 15
points = 101
x = np.linspace(-v, v, points)
y = np.linspace(-v, v, points)
pix_size = (2*v)/points

#have to replace this with phenomenological relations
rt = [i*2 for i in rh]
rs = rvmax

t = tau(rs,rt)
mnfw = mnfw(mass2,t)

count = 0
conv_list = []
for g in pos_2d:
    count += 1
    print count
    result = [conv_tnfw(x[:,None]+i[0],y[None,:]+i[1],k,l,m).real for i,k,l,m in zip(g,t,mnfw,rs) if -30 < i[0] < 30 and -30 < i[1] < 30]
    res = np.log10(sum(result))
    conv_list.append(res)

conv_array = []
for i in conv_list:
    k = [np.asarray(j) for j in i]
    conv_array.append(k)

avg_conv = [sum(f)/len(conv_array) for f in zip(*conv_array)]

conv_list2 = np.asarray(conv_list)
avg_conv2 = np.asarray(avg_conv)

np.save('kavg',np.mean(kavg))
np.save('ind_conv',conv_list2)
np.save('tot_conv',avg_conv2)

print("--- %s seconds ---" % (time.time() - start_time))

from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
from R_ein import rein_sigmac
from functions import simul_data,rotation,projections,angular_average,twoD_ps
from scipy.optimize import curve_fit

start_time = time.time()

posmass = np.load('projections.npy')

rnge = 100
rein = 10.300324225
sigmac = 3120194565.15

pos_2d = []
masses = []
for i in posmass:
    r,m,_,_ = zip(*i)
    pos_2d.append(r)
    masses.append(m)    

num = [len(i) for i in masses]
minmass = [min(i) for i in masses]
maxmass = [max(i) for i in masses]
mean_mass = [np.mean(i) for i in masses]

print 'average num of subh = %s ' % np.mean(num)
print 'average msub = %s' % np.mean(mean_mass)
print 'average min subhalo mass = %s +- %s' % (np.mean(minmass),np.std(minmass))
print 'average max subhalo mass = %s +- %s' % (np.mean(maxmass),np.std(maxmass))

flat_list = [item for sublist in masses for item in sublist]

N,bins,bars = plt.hist(flat_list,log=True,bins=np.logspace(6.0,8.0,75))

dm = []
count = 0
for i in range(len(bins)-1):
    count += 1
    dm.append(bins[count] - bins[count-1])

tab = [[i/len(posmass),j,k] for i,j,k in zip(N,bins,dm) if i != 0]
N,bins,dm = zip(*tab)

dNdm = [i/j for i,j in zip(N,dm)]

Nsub = np.sum([i*j for i,j in zip(dm,dNdm)])
print 'recovering the average number of subhalos from the mass function: %s' % Nsub

def func(x, a0):
    return a0 * (x/(2.52*10**7))**(-1.9)

params, pcov = curve_fit(func, bins, dNdm)

print 'best fit params %s' % params
f = [func(i,params[0]) for i in bins]


"""#a0 = 6.98854*1e-6
a0 = 2.577227*1e-7
m0 = 2.52*1e7
beta = -1.9
f = [func(i,a0) for i in bins]"""


#aconv = 2.6488187356180315*1e-6
#f3 = [func(i,aconv) for i in bins]

plt.rc('text', usetex=True)

plt.loglog(bins,dNdm,label='average $\\frac{dN}{dm}$ from projections')
#plt.loglog(bins,f,label='a0 = %s ' % (a0),c='b')
plt.loglog(bins,f,label='a0 = %s' % (params[0]),c='g')
#plt.loglog(bins,f3,label='a0 = %s' % (aconv),c='m')
plt.ylabel('dN/dM [$M_{\odot}^{-1}$]')
plt.xlabel('Mass [$M_{\odot}$]')
plt.xlim(min(bins),max(bins))
plt.ylim(min(dNdm),max(dNdm))
plt.legend()
plt.show()

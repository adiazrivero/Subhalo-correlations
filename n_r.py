from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
from functions import simul_data,rotation,projections,poisson_realization

start_time = time.time()

list1 = []
for i in range(16):
    list1.append('/Users/anadiazrivero/Desktop/ETHOS/CDM/fof_subhalo_tab_127.%s.hdf5' % i)

Rng = 600
subh_mass,subh_pos,r1,r2 = simul_data(list1,Rnge=Rng) #in 10^10 M_sun/h, kpc/h
subh_mass = [i*10**(10) for i in subh_mass]

plot_N3d = True
plot_N2d = False

if plot_N3d == True:
    r = [np.sqrt(i**2+j**2+k**2) for i,j,k, in subh_pos]
    binsr = np.arange(1,500,5)
    counts,bins,bars = plt.hist(r,bins=binsr)
    plt.ylabel('N_sub')
    plt.xlabel('r_3D (kpc/h)')
    plt.show()
    plt.clf()

def nr_3d(coords,masses,mcut=False,mhigh=1,mlow=1):

    if mcut == True:
        arr = zip(coords,masses)
        arr2 = []
        for a,b in zip(mlow,mhigh):
            arr3 = [i for i in arr if a < i[1] < b]
            arr2.append(arr3)
        x = []
        y = []
        z = []
        for i in arr2:
            pos,mass = zip(*i)
            xx,yy,zz = zip(*pos)
            x.append(xx)
            y.append(yy)
            z.append(zz)

        #print np.shape(x),np.shape(y),np.shape(z)

    else:
        x,y,z = zip(*coords)

    pix_num = 50
    pix_size = Rng/pix_num
    H = []
    xe = []
    ye = []
    ze = []
    Rmax = []
    for i,j,k in zip(x,y,z):
        i = np.asarray(i)
        j = np.asarray(j)
        k = np.asarray(k)
        r2 = np.sqrt(i**2+j**2+k**2)
        Rmax.append(max(r2))
        counts,[xedges,yedges,zedges] = np.histogramdd((i,j,k),bins=(pix_num,pix_num,pix_num))
        H.append(counts)
        xe.append(xedges)
        ye.append(yedges)
        ze.append(zedges)

    xe2 = []
    for i in xe:
        xedges = [j+0.5*pix_size for j in i[:-1]]
        xe2.append(xedges)

    ye2 = []
    for i in ye:
        yedges = [j+0.5*pix_size for j in i[:-1]]
        ye2.append(yedges)

    ze2 = []
    for i in ze:
        zedges = [j+0.5*pix_size for j in i[:-1]]
        ze2.append(zedges)

    nr_tot = []
    R_tot = []
    for i,j,k,l,q in zip(H,xe2,ye2,ze2,Rmax):
        im = np.array(i)
        xx,yy,zz = np.meshgrid(j,k,l)
        r = np.sqrt(xx**2+yy**2+zz**2)
        rmax = q
        dr = pix_size
        R = np.arange(rmax/dr)*dr
        R_tot.append(R)
        nr = []
        for u in range(len(R)):
            rmin = u*dr
            rmax = rmin + dr
            index = (r>=rmin) * (r<=rmax)
            nr.append(im[index].mean())
        nr_tot.append(nr)

    return R_tot,nr_tot

mhigh = [1e10,1e9,1e8,1e7]
mlow = [1e9,1e8,1e7,1e6]
R,nr = nr_3d(subh_pos,subh_mass,mcut=True,mhigh=mhigh,mlow=mlow)

plt.rc('text', usetex=True)

mhigh = [10,9,8,7]
mlow = [9,8,7,6]
for i,j,k,l in zip(nr,R,mlow,mhigh):
    aq = np.asarray(i)/np.mean(i)
    #plt.loglog(j,aq,label='$%.1e M_{\odot} < M_{sub} < %.1e M_{\odot}$' % (k,l))
    plt.loglog(j,aq,label='$10^{%s} M_{\odot} < M_{sub} < 10^{%s} M_{\odot}$' % (k,l))
    plt.xlabel('$r_{3D}$ $[kpc/h]$',fontsize=20)
    plt.ylabel('$n(r)/<n>$',fontsize=20)

#plt.title('n(r) in 3d')
plt.tick_params(axis='both', which='major', labelsize=15)
plt.legend(prop={'size': 16})
plt.show()

shft = 0
rnge = 1000
coordsm,avg_num_subh = projections(subh_pos,subh_mass,r1,r2,rnge=rnge,shift=shft,num_proj=100)
#coordsm,avg_num_subh = projections(subh_pos,subh_mass,r1,r2,rnge=rnge,shift=shft,num_proj=100)
print "Done"

coords = []
for i in coordsm:
    r,m,_,_ = zip(*i)
    coords.append(r)

if plot_N2d == True:
    #plot as 1d hist
    pos_2d = [i for j in coords for i in j]
    r_2d = [np.sqrt(i**2+j**2) for i,j in pos_2d]
    binsr = np.arange(1,600,5)
    counts,bins,bars = plt.hist(r_2d,bins=binsr)
    plt.ylabel('N_sub')
    plt.xlabel('r_2D (kpc/h)')
    plt.show()
    plt.clf()
    #plot as line
    counts = [i/len(coords) for i in counts]
    diff = binsr[1]-binsr[0]
    binsr = [i + 0.5*diff for i in binsr][:-1]
    plt.plot(binsr,counts)
    plt.ylabel('N_sub')
    plt.xlabel('r_2D (kpc/h)')
    plt.show()

def nr_2d(coords,bns=20,rnge=100):
    bin_size = rnge/bns
    nr = []
    n_avg = []
    xlim = []
    ylim = []
    for i in coords:
        x,y = zip(*i)
        N,xedges,yedges = np.histogram2d(x,y,bins=bns)
        xedges = [i+0.5*bin_size for i in xedges[:-1]]
        yedges = [i+0.5*bin_size for i in yedges[:-1]]
        xlim.append(xedges)
        ylim.append(yedges)
        n = [i/bin_size**2 for i in N]
        nr.append(n)
        Nbar = N.mean()
        n_avg.append(Nbar/bin_size**2)
    return nr,n_avg,bin_size,xlim,ylim

nr,n_avg,bin_size,xlim,ylim = nr_2d(coords,bns=500,rnge=rnge)

print bin_size

u = [sum(x)/len(coords) for x in zip(*nr)]
xavg = [sum(i)/len(coords) for i in zip(*xlim)]
yavg = [sum(i)/len(coords) for i in zip(*ylim)]

im = np.array(u)
xx,yy = np.meshgrid(xavg,yavg)
r = np.sqrt(xx**2+yy**2)
rmax = max(xavg)
dr = bin_size
R = np.arange(rmax/dr)*dr

nr = []
for i in range(len(R)):
    rmin = i*dr
    rmax = rmin + dr
    index = (r>=rmin) * (r<=rmax)
    nr.append(im[index].mean())

plt.rc('text', usetex=True)

fig,ax = plt.subplots(1)
ax.axhline(nr[1],min(R),1000,c='r')
ax.axvline(15,1e-5,1,c='r',ls='dashed')
ax.loglog(R,nr)
ax.set_xlabel('$r_{2D}$ $[kpc/h]$',fontsize=20)
ax.set_ylabel('$n(r)$',fontsize=20)
ax.set_xlim(min(R),1000)
ax.set_ylim(min(nr),2*1e-1)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.title('n(r) in 2d')
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

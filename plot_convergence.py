from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pylab as py
import sys
import time
import cmath
from rein import rein_sigmac
from functions import simul_data,rotation,projections,twoD_ps,angular_average

kavg_14 = np.load('kavg.npy')
conv_list = np.load('ind_conv.npy')
avg_conv = np.load('tot_conv.npy')

print kavg_14

rein,sigmac = rein_sigmac(1.04e12,1,0.5) #in kpc and M_sun/kpc^2
h = 0.7
rein_kpch = h*rein
num_proj = 30

rnge = 100
shft = 0
v = 15
points = 101
x = np.linspace(-v, v, points)
y = np.linspace(-v, v, points)
pix_size = (2*v)/points

###################################################################################
#finding average convergence at r_ein and the average convergence in the SL region
###################################################################################

im = np.array(avg_conv)
im2 = 10**im
xx,yy = np.meshgrid(x,y)
r = np.sqrt(xx**2+yy**2)
rmax = v
dr = pix_size
R = np.arange(rmax/dr)*dr

kappar = []
kappa_rein = []
cum_avgkappa = []

for i in range(len(R)):
    rmin = i*dr
    rmax = rmin + dr
    index = (r>=rmin) * (r<=rmax)
    if rein-0.5 < rmax < rein+0.5:
        kappa_rein.append(im2[index].mean())
    if rmax <= rein:
        cum_avgkappa.append(im2[index].mean())

avgkappa = np.mean(kappa_rein)


print "Average convergence at the Einstein radius: %s " % avgkappa
print "Average convergence in the SL region: %s " % np.mean(cum_avgkappa)

###################################################################################
#plot the 2d convergence maps
###################################################################################

plot_2dkappa = False

if plot_2dkappa == True:
    fig,ax = plt.subplots(1)
    k = ax.imshow(avg_conv,extent=[-v,v,-v,v])
    Rein = rein
    circ = Circle((0,0),Rein,facecolor='none',fill=False,linestyle='dashed')
    ax.add_patch(circ)
    fig.colorbar(k)
    plt.xlabel('kpc/h')
    plt.title('Convergence field on the lens plane (%s proj/kavg=%.4f/kavg_14 = %.4f)' % (num_proj*3,np.mean(cum_avgkappa),kavg_14))
    plt.show()

####################################################################
# Doing the 2D FT and angular-averaging to obtain the PS
####################################################################

individual_ps,tot_ps,kx,ky = twoD_ps(data=None,ind_data=conv_list,pix_size=pix_size,rnge=rnge,shift=shft,show_ps=False)
pix_size_k = np.abs(kx[0]-kx[1])
ps1d,K = angular_average(tot_ps,kx,ky,rnge=rnge,pix_num=points,dr=pix_size_k,remove_first=True)

py.figure('2')
plt.title('Power Spectrum of the convergence field (%s proj/kavg=%.4f/kavg_14=%.4f)' % (num_proj*3,np.mean(cum_avgkappa),kavg_14))
plt.loglog(K,ps1d,label='ind')
plt.xlabel('k [h kpc]^-1')
plt.ylabel('P(k) [kpc/h]^2')
plt.legend()
plt.show()
plt.clf()

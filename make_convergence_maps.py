from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pylab as py
import sys
import time
#from rein import rein_sigmac
from functions import simul_data,rotation,projections,angular_average,twoD_ps


posmass = np.load('projections.npy')
m_high_cut = True

rein= 10.300324225
sigmac = 3120194565.15

pos_2d = []
masses = []
for i in posmass:
    r,m,_,_ = zip(*i)
    pos_2d.append(r)
    masses.append(m)

def rs(msub):
    return 0.1*(msub/10**6)**(1/3)

def mnfw(m,tau):
    n = np.asarray([(i**2/(i**2 + 1)**2) * ((i**2 - 1) * np.log(i) + i * np.pi - (i**2 + 1)) for i in tau])
    mnfw = m/n
    return mnfw

def FNFW(x):
    up = x>1
    down = np.logical_not(up)
    lowx = np.arccosh(1/x[down])/np.sqrt(1 - x[down]**2)
    highx = np.arccos(1/x[up])/np.sqrt(x[up]**2 - 1)
    ans = np.zeros(x.shape)
    ans[down] = lowx
    ans[up] = highx
    return ans

def conv_tnfw(x,y,tau,mnfw,rs):
    z = np.sqrt(x**2+y**2)/rs
    Lz = np.log(z/(tau+np.sqrt(tau**2+z**2)))
    Fz = FNFW(z)
    sig = (mnfw/(sigmac*rs**2))*(tau**2/(2*np.pi*(tau**2 + 1)**2))*(((tau**2 + 1)/(z**2 - 1))*(1 - Fz) + 2*Fz - np.pi/np.sqrt(tau**2 + z**2) + ((tau**2 - 1)/(tau*np.sqrt(tau**2 + z**2)))*Lz)
    #print(np.argwhere(np.isnan(sig)))
    return sig

masses = np.asarray([np.asarray(i) for i in masses])
r_s = rs(masses)
avg_min_rs = np.mean([min(i) for i in r_s])

rt = 5*r_s
max_rt = [max(i) for i in rt]
avg_max_rt = np.mean(max_rt)

t = rt/r_s
m_nfw = mnfw(masses,t)

rnge = 100
shiftx = 0
shifty = 0
pix_num = 101
pix_size = rnge/pix_num
x = np.linspace(-rnge/2+shiftx, rnge/2+shiftx, pix_num)
y = np.linspace(-rnge/2+shifty, rnge/2+shifty, pix_num)

count = 0
conv_list = []

for g,q,u,z in zip(pos_2d,t,m_nfw,r_s):
    count += 1
    print count
    beginning = time.time()
    result = [conv_tnfw(x[:,None]-i[0],y[None,:]-i[1],k,l,m).real for i,k,l,m in zip(g,q,u,z)]
    res = sum(result)
    conv_list.append(res)
    print("%s seconds" % (time.time() - beginning))

    if count % 10 == 0:
        conv_list2 = np.asarray(conv_list)
        np.save('conv_%s_%s_%s' % (pix_num,count,rnge),conv_list2)
        conv_list = []

feat = {'pixels':pix_num,'rnge':rnge,'num_proj':len(posmass),'m_high cut':m_high_cut}

np.save('convfeat_%s_%s' % (pix_num,rnge), feat)

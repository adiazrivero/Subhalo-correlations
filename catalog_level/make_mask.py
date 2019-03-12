from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
from functions import simul_data,rotation,projections,twoD_nr,twoD_ps,variance,poisson_realization,interpolation

def make_mask(kx,ky,dk=1):
    x,y = np.meshgrid(ky,kx)
    k = np.sqrt(x**2+y**2)
    kmax = max(kx)
    dK = dk
    K = np.arange(kmax/dK)*dK
    mask_list =[]
    for i in range(len(K)):
        kmin = i*dK
        kmax = kmin + dK
        mask = (k >= kmin) * (k <= kmax)
        mask_list.append(mask*1)
        #print np.count_nonzero(mask*1)
        #print mask*1
    return mask_list

pix_num = 1011
rnge = 100
pix_size = rnge/pix_num
X = np.linspace(-rnge/2.,rnge/2.,pix_num)
Y = np.linspace(-rnge/2.,rnge/2.,pix_num)
mask_list = make_mask(X,Y,dk=pix_size)

#print np.sum(mask_list,axis=0)
np.save('mask_%s.npy' % pix_num, mask_list)

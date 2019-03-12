from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import argparse
#from functions import simul_data,rotation,projections,twoD_nr,twoD_ps,variance,poisson_realization,interpolation

parser = argparse.ArgumentParser()
parser.add_argument("-p","--pix_num",
                    default=1011,
                    help="number of pixels in the image",
                    type=int)
parser.add_argument("-s","--side",
                    default=100,
                    help='physical size of the image in kpc',
                    type=int)

args = parser.parse_args()

pix_num = args.pix_num
rnge = args.side

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
    return mask_list

#rnge = 100
pix_size = rnge/pix_num
X = np.linspace(-rnge/2.,rnge/2.,pix_num)
Y = np.linspace(-rnge/2.,rnge/2.,pix_num)
mask_list = make_mask(X,Y,dk=pix_size)

np.save('mask_%s.npy' % pix_num, mask_list)

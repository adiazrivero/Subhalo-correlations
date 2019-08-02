from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import cmath
import argparse
from functions import *

h = 0.6909

parser = argparse.ArgumentParser()

parser.add_argument('--conv_file1',
                  help='path for file with convergence maps, output of convergence_maps.py',
                  type=str)
parser.add_argument('--conv_file2',
                  help='path for file with convergence maps, output of convergence_maps.py',
                  type=str)
parser.add_argument('--conv_file3',
                  help='path for file with convergence maps, output of convergence_maps.py',
                  type=str)
parser.add_argument('--kdir',
                    default='./',
                    help='output directory for wavenumbers',
                    type=str)
parser.add_argument('--psdir',
                    default='./',
                    help='output directory for the power spectra',
                    type=str)
parser.add_argument("-p","--pix_num",
                    default=1011,
                    help="number of pixels in the image",
                    type=int)
parser.add_argument("-s","--side",
                    default=100,
                    help='physical size of the image in kpc',
                    type=int)
parser.add_argument("--name",
                    default='CDM',
                    help='which DM model to use',
                    type=str)
parser.add_argument("-z","--z",
                    default=0.5,
                    help="simulation redshift",
                    type=float)
parser.add_argument("--m_high",
                   default=1e-2,
                   help='highest subhalo mass, in units of 10^10 M_sun',
                   type=float)
parser.add_argument("--m_low",
                   default=1e-4,
                   help='highest subhalo mass, in units of 10^10 M_sun',
                   type=float)
#parser.add_argument("-n","--num_proj",
                   #default=10,
                   #help='total number of projections divided by 3',
                   #type=int)
#parser.add_argument('-S','--sigma_crit',
                    #default=2.35e9,
                    #help='critical surface mass density for lensing, in units of M_sun/kpc^2',
                    #type=float)

args = parser.parse_args()
conv_file1 = args.conv_file1
conv_file2 = args.conv_file2
conv_file3 = args.conv_file3
dirk = args.kdir
dirps = args.psdir
pix_num = args.pix_num
rnge = args.side
name = args.name
z = args.z
mhigh = args.m_high
mlow = args.m_low
#num_proj = args.num_proj
#sigmac = args.sigma_crit

######################################################

if z == 0.5:
    numb = '095'

elif z == 0:
    numb='127'
    z = 0.5 #redshift changed to give physical values  

elif z == 1:
    numb = '078'

else:
    print("wrong redshift specified")
    sys.exit()

######################################################

if mhigh >= 10:
    mlab = 'm11'

elif mhigh == 1:
    mlab = 'm10'

elif mhigh ==1e-1:
    mlab = 'm9'

elif mhigh == 1e-2:
     mlab = 'm8'

elif mhigh == 1e-3:
    mlab = 'm7'

######################################################


if mlow == 1e-4:

    mpart = 1.90382332e-6 / h
    m_min = 50 * mpart
    mlab2 = 'm6'

elif mlow == 1e-3:

    mlab2 = 'm7'
    m_min = mlow

elif mlow == 1e-2:

    mlab2 = 'm8'
    m_min = mlow

elif mlow == 1e-1:

    mlab2 = 'm9'
    m_min = mlow

elif mlow == 1:

    mlab2 = 'm10'
    m_min = mlow

elif mlow == None:
    mpart = 1.90382332e4 / h
    m_min = 50 * mpart
    mlab2 = 'm6'
    m_min = mlow

if mlab == mlab2 or mhigh <= m_min:
    print('wrong mass bounds specified!')
    sys.exit()

######################################################
start = time.time()
############################################
#loading the convergence maps
############################################

pix_size = rnge/pix_num

conv_list1 = np.load(conv_file1)
conv_list2 = np.load(conv_file2)
conv_list3 = np.load(conv_file3)

print('convergence maps loaded:')
print(conv_file1)
print(conv_file2)
print(conv_file3)

conv_list = np.concatenate((conv_list1,conv_list2,conv_list3))

####################################################################
# Doing the 2D FT and angular-averaging to obtain the PS
####################################################################

individual_ps,tot_ps,kx,ky = twoD_ps(data=conv_list,pix_size=pix_size,rnge=rnge)

#load a mask and keep only the right number of rings (depends on pixel number)

orig_mask = np.load('mask_1011.npy')
bin_num = int((pix_num)/2 + 1)
mask = []
for i in orig_mask[:bin_num]:
    mask2 = []
    beg = int((len(i[0])-pix_num)/2)
    end = int(len(i[0])-beg)
    for j in i[beg:end]:
        mask2.append(j[beg:end])
    mask.append(mask2)

pix_size_k = np.abs(kx[0]-kx[1])

#obtain the 1d power spectrum for each convergence map
ns = [0]
power_spectra,K = multipoles(tot_ps,kx,ky,mask=mask,pix_num=pix_num,dr=pix_size_k,ns=ns)

K = K[1:]

dirk = dirk + 'k_%s_%s_%s_%s.txt' % (name,numb,pix_num,rnge)
file0 = open(dirk,'w')
for i in K:
    file0.write('%s\n' % i)
file0.close()
print('k saved in %s' % dirk)

for key in power_spectra.keys():
    
    ps = power_spectra[key][1:] 
    _,ind_curves = variance(individual_ps,ps,len(individual_ps),kx,ky,mask=mask,rnge=rnge,pix_num=pix_num,pix_size=pix_size_k,n=key)

    dir2 = dirps + 'ind_curves_%s_%s_%s_%s' % (name,numb,pix_num,rnge)
    np.save(dir2,ind_curves)

    print('saved %s' % dir2)

end = time.time() - start
print(end)

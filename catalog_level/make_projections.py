from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from functions import *
import cPickle as pickle
import argparse

start_time = time.time()

h = 0.6909

parser = argparse.ArgumentParser()

parser.add_argument('-o','--outdir',
                    default='./',
                    help='output directory',
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
                    help="simulation redshift")
parser.add_argument("--m_high",
                   default=1e-2,
                   help='highest subhalo mass')
parser.add_argument("--m_low",
                   default=1e-4,
                   help='highest subhalo mass')
parser.add_argument("-n","--num_proj",
                   default=10,
                   help='total number of projections divided by 3')

args = parser.parse_args()

outdir = args.outdir
pix_num = args.pix_num
rnge = args.side
name = args.name
z = args.z
mhigh = args.m_high
mlow = args.m_low
num_proj = args.num_proj

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

#####################################################

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

sys.exit()

######################################################

list1 = []
for i in range(16):
    list1.append('/n/hernquistfs3/jzavala/ETHOS/%s/groups_%s/fof_subhalo_tab_%s.%s.hdf5' % (name,numb,numb,i))    

shft = 0

#obtaining subhalo catalogs

m_min = m_min * 1e10
m_high = mhigh * 1e10
print 'lower mass bound: %.3e' % m_min
print 'upper mass bound: %.3e' % m_high

mass,pos,rh,rvmax = simul_data(list1,mlow=m_min,mhigh=m_high,Rnge=600)

print "total number of subhalos within L_box =  300 kpc: %s " % np.shape(mass)
print "max mass = %e" % max(mass)
print "min mass = %e" %  min(mass)

posmass,_ = projections(pos,mass,rh,rvmax,rnge=rnge,shift=shft,num_proj=num_proj)
posmass2 = np.asarray(posmass)
print "proj num = %s " % len(posmass2)

pos_2d = []
for i in posmass:
    r,_,_,_,_ = zip(*i)
    pos_2d.append(r)

N_sub = [len(i) for i in pos_2d]

med = np.median(N_sub)
plus = np.percentile(N_sub,95) - med
minus = med - np.percentile(N_sub,5)
print 'N_sub = %.2f + %.2f - %.2f' % (med,plus,minus)

filename = outdir + 'projections_%s_%s_%s' % (name,numb,rnge)

np.save(filename,posmass2)
print 'saved as %s' % filename


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from functions import *
import cPickle as pickle

start_time = time.time()

#executing variables declared externally

list=sys.argv[1:]
for i in list:
    exec(i)

if numb == 78 or numb == 95:
    numb = '0' + str(numb)
else:
    numb = str(numb)

if name == 0:
    name = 'CDM'
elif name == 1:
    name = 'ETHOS_1'
elif name == 2:
    name = 'ETHOS_4'

print name

if numb == '127':
    z = 0.5
    print "z = %s; redshift changed to give physical values" % z
elif numb == '095':
    z = 0.5
    print 'z = %s' % z
elif numb == '078':
    z = 1
    print 'z = %s' % z
else:
    print "wrong redshift specified"
    sys.exit()

list1 = []
for i in range(16):
    list1.append('/n/hernquistfs3/jzavala/ETHOS/%s/groups_%s/fof_subhalo_tab_%s.%s.hdf5' % (name,numb,numb,i))    

shft = 0
rnge = 100
num_proj = 10 #actual number of projections is this * 3

#obtaining subhalo catalogs

mass,pos,rh,rvmax = simul_data(list1,mhigh_cut=m_high_cut,mhigh=mhigh,Rnge=600) #mhigh is in units of 10^10 M_sun

print "total number of subhalos within r_3d < 300 kpc: %s " % np.shape(mass)

mass2 = [i*10**(10) for i in mass] #units of M_sun

print "max mass = %e" % max(mass2)
print "min mass = %e" %  min(mass2)

dicti = pickle.load(open('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/host_features/host_%s_%s.txt' % (name,numb),'rb'))

Mhost = dicti['host_mass']

arr = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/rein_sigmac_all/rein_sigmac_%s_%s.npy' % (name,numb))
zs = arr[0]
zl = arr[1]
rein = arr[3]
sigmac = arr[4]

if zl != z:
    print "The wrong simulation is being used as a lens!"
    sys.exit()

print "r_ein = %s kpc " % rein
print "sigma_crit = %e M_sun/kpc^2 " % sigmac

posmass,_ = projections(pos,mass2,rh,rvmax,rnge=rnge,shift=shft,num_proj=num_proj)

posmass2 = np.asarray(posmass)
print "proj num = %s " % len(posmass2)

if mhigh == 0:
    np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/m0/projections_%s_%s_m0' % (name,numb,name,numb), posmass2)
elif mhigh == 1e-2:
    np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/m8/projections_%s_%s_m8' % (name,numb,name,numb), posmass2)
elif mhigh == 1e-3:
    np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/m7/projections_%s_%s_m7' % (name,numb,name,numb), posmass2)


from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from nbodykit.source.catalog import HDFCatalog,CSVCatalog
from nbodykit import transform
from nbodykit.lab import ProjectedFFTPower, FFTPower,BigFileMesh
from nbodykit.source.mesh.array import ArrayMesh
import dask.array as da
import cPickle as pickle
from tnfw_convergence import *
from tburk_convergence import *
sys.path.insert(0, '/n/home04/adiazrivero/conv_powerspec/catalog_level/')
from functions import *

bs = 600
nmesh = 1024
h = 0.6909
mpart = 1.90382332e4 / h

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

arr = np.load('/n/home04/adiazrivero/conv_powerspec/catalog_level/rein_sigmac_all/rein_sigmac_%s_%s.npy' % (name,numb))
zs = arr[0]
zl = arr[1]
rein = arr[3]
sigmac = arr[4]

#cosmo_lens = LensCosmo(zl,zs)
#sigmac = cosmo_lens.sigma_crit/(cosmo_lens.arcsec_2_kpc)**2
print zl,zs,'%.4e' % sigmac
#print "r_ein = %.2f kpc " % rein
#print "sigma_crit = %.2e M_sun/kpc^2 " % sigmac

if zl != z:
    print "The wrong simulation is being used as a lens!"
    sys.exit()

dicti = pickle.load(open('/n/home04/adiazrivero/conv_powerspec/particle_level/nbar/hbugfix/nbar_%s_%s.txt' % (name,numb),'rb'))
N_part = dicti['N_part']
nbar = dicti['nbar'] #for 600 kpc box

print 'nbar = %.3f' % nbar

#mesh = BigFileMesh('/n/dvorkin_lab/anadr/%s_%s_%s_%s_h0.6909.bigfile' % (name,numb,nmesh,bs),'Field')
meshdir = '/n/dvorkin_lab/anadr/%s_%s_%s_%s_hbugfix.bigfile' % (name,numb,nmesh,bs)
mesh = BigFileMesh(meshdir,'Field')

print 'loaded mesh from %s' % meshdir

one_plus_delta = mesh.paint(mode='real') #compensate for wrong h before
l = bs/nmesh
V = l**3
n = one_plus_delta.value * nbar 
N = n * V

print('input number of particles:%s // recovered number: %s' % (N_part,N.sum()))

##################################################
#average map \approx host
##################################################

rho = n * mpart

m_tot = np.sum(rho) * l**2
print 'total mass in projection: %.4e' % m_tot

rho = ArrayMesh(rho,(bs,bs,bs),Nmesh=nmesh)
kappa0 = (1/sigmac) * rho.preview(axes=[0,1], Nmesh=nmesh)

print 'kappa0 max = %s' % np.amax(kappa0)

np.save('kappa0_%s_%s_%s_%s' % (bs,name,numb,nmesh),kappa0)

kappa0 = np.fft.fftshift(kappa0) #[426:-426,426:-426] #trimming box so that it is only +- 50 kpc 
#bs = 100

print 'kappa0 shape = %s' % kappa0.shape[0]
print 'kappa0 max = %s' % np.amax(kappa0)

##################################################
#single projection without subtracting host
##################################################

kappa0_mesh = ArrayMesh(kappa0,(bs,bs),Nmesh=kappa0.shape)

Pk0 = ProjectedFFTPower(kappa0_mesh,axes=(0,1),BoxSize=bs,Nmesh=kappa0.shape)
ks0 = Pk0.power['k']
ps0 = Pk0.power['power']

file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/singleprojp_%s_%s_%s.txt'% (name,numb,nmesh), 'w')
for i in ps0.real:
    file1.write('%s\n' % i)
file1.close()
file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/singleprojk_%s_%s_%s.txt' % (name,numb,nmesh), 'w')
for i in ks0:
    file1.write('%s\n' % i)
file1.close()

##################################################
#estimating the host
##################################################

allprojs = []
allprojs.append(kappa0)

"""for i in range(1,201):
    kappa1 =  nbar * mpart / sigmac * np.load('/n/dvorkin_lab/anadr/%s_proj/%s/%s_%s_projxy_%s.npy' % (name,numb,name,numb,i)) / h**2  #[426:-426,426:-426]
    kappa2 =  nbar * mpart / sigmac * np.load('/n/dvorkin_lab/anadr/%s_proj/%s/%s_%s_projxz_%s.npy' % (name,numb,name,numb,i)) / h**2 #[426:-426,426:-426]
    kappa3 =  nbar * mpart / sigmac * np.load('/n/dvorkin_lab/anadr/%s_proj/%s/%s_%s_projyz_%s.npy' % (name,numb,name,numb,i)) / h**2 #[426:-426,426:-426]
    
    #print 'kappa1 max = %s' % np.amax(kappa1)
    
    allprojs.append(kappa1)
    allprojs.append(kappa2)
    allprojs.append(kappa3)"""

for i in range(1,11):
    for j in range(1,21):
        kappa1 =  nbar * mpart / sigmac * np.load('/n/dvorkin_lab/anadr/%s_proj/%s/hbugfix/%s_%s_projxy_%s_%s.npy' % (name,numb,name,numb,i,j))  #[426:-426,426:-426]
        kappa2 =  nbar * mpart / sigmac * np.load('/n/dvorkin_lab/anadr/%s_proj/%s/hbugfix/%s_%s_projxz_%s_%s.npy' % (name,numb,name,numb,i,j))  #[426:-426,426:-426]
        kappa3 =  nbar * mpart / sigmac * np.load('/n/dvorkin_lab/anadr/%s_proj/%s/hbugfix/%s_%s_projyz_%s_%s.npy' % (name,numb,name,numb,i,j))  #[426:-426,426:-426]
    
        #print 'kappa1 max = %s' % np.amax(kappa1)
    
        allprojs.append(kappa1)
        allprojs.append(kappa2)
        allprojs.append(kappa3)

print len(allprojs)

print 'kappa1 shape = %s' % kappa1.shape[0]
print 'kappa1 max = %s' % np.amax(kappa1)

"""kappa1_mesh = ArrayMesh(kappa1,(bs,bs),Nmesh=kappa0.shape[0])
Pk0 = ProjectedFFTPower(kappa1_mesh,axes=(0,1),BoxSize=bs,Nmesh=kappa0.shape[0])
ks0 = Pk0.power['k']
ps0 = Pk0.power['power']

file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/singleprojp_%s_%s.txt'% (name,numb), 'w')
for i in ps0.real:
    file1.write('%s\n' % i)
file1.close()
file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/singleprojk_%s_%s.txt' % (name,numb), 'w')
for i in ks0:
    file1.write('%s\n' % i)
file1.close()

sys.exit()"""

host = np.mean(allprojs,axis=0)

print 'host shape = %s' % host.shape[0]

np.save('host_%s_%s_%s' % (bs,name,numb),host)

print 'host amplitude: %s' % np.amax(host)

host_mesh = ArrayMesh(host,(bs,bs),Nmesh=host.shape[0])
Pk1 = ProjectedFFTPower(host_mesh,axes=(0,1),BoxSize=bs,Nmesh=host.shape[0])
#plt.loglog(Pk1.power['k'],Pk1.power['power'].real*ratio,label='average without subtracting host')

file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/avgwhostp_%s_%s.txt' % (name,numb), 'w')
for i in Pk1.power['power'].real:
    file1.write('%s\n' % i)
file1.close()
file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/avgwhostk_%s_%s.txt' % (name,numb), 'w')
for i in Pk1.power['k']:
    file1.write('%s\n' % i)
file1.close()

#################################################################
#single projection subtracting the host
#################################################################

sub = kappa0 - host
sub1 = kappa1 - host 
sub2 = kappa2- host 
sub3 = kappa3 - host 

np.save('sub1_%s_%s_%s' % (bs,name,numb),sub)
np.save('sub2_%s_%s_%s' % (bs,name,numb),sub2)
np.save('sub3_%s_%s_%s' % (bs,name,numb),sub3)

sub_mesh = ArrayMesh(sub,(bs,bs),Nmesh=sub.shape[0])
Pksub = ProjectedFFTPower(sub_mesh,axes=(0,1),BoxSize=bs,Nmesh=sub.shape[0])

ksub = Pksub.power['k']
pksub = Pksub.power['power']

#plt.loglog(ksub,pksub.real*ratio,label='single subtracting host')

file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/singlewohostp_%s_%s.txt'% (name,numb), 'w')
for i in pksub.real:
    file1.write('%s\n' % i)
file1.close()
file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/singlewohostk_%s_%s.txt' % (name,numb), 'w')
for i in ksub:
    file1.write('%s\n' % i)
file1.close()

#################################################################
# fitting tNFW profile to host
#################################################################

"""pix_size = bs/nmesh
x = np.linspace(-bs/2, bs/2, nmesh)
y = np.linspace(-bs/2, bs/2, nmesh)

dicti = pickle.load(open('/n/home04/adiazrivero/conv_powerspec/catalog_level/host_features/host_%s_%s.txt' % (name,numb),'rb'))
Mhost = dicti['host_mass']
Rmaxhost = dicti['host_rmax']
print '%e' % Mhost
print '%s' % Rmaxhost
sys.exit()

if name == 'CDM':

    fit_lab = 'tnfw'
    #r_s = 63/2.1626
    r_s = 10
    r_200 = 244.05
    t = r_200/r_s
    Mhost = 1.6128e12
    m_nfw = mnfw(Mhost,t)

    print 'm_nfw = %e' % m_nfw

    fit = conv_tnfw(x[:,None],y[None,:],t,m_nfw,r_s,sigmac).real

    print 'fit amplitude: %s' % np.amax(fit)

elif name == 'ETHOS_4':
    
    fit_lab = 'tburk'
    #r_s = 69.18/2.1626
    r_s = 20
    r_200 = 245.3
    p = 0.2
    #p = 0.0625
    r_core = p * r_s
    print 'core radius = %.3f' % r_core
 
    t = r_200/r_s
    Mhost = 1.6376e12
    m_b = mb(Mhost,t,p)

    print 'm_b = %e' % m_b

    fit = conv_tburk_2d(x[:,None],y[None,:],r_s,m_b,p,t,sigmac).real

    print 'fit amplitude: %s' % np.amax(fit)


smd = fit*sigmac

m_tot = np.sum(smd)

print 'total mass in projection = %.3e' % m_tot

np.save('%s_fit_convergence_%s_%s_2' % (fit_lab,name,numb),fit)
 
#################################################################
#single projection subtracting the fit profile
#################################################################

sub = kappa1 - fit
sub_mesh = ArrayMesh(sub,(bs,bs))
Pksub = ProjectedFFTPower(sub_mesh,axes=(0,1),BoxSize=bs,Nmesh=nmesh)

ksub = Pksub.power['k']
pksub = Pksub.power['power']

#plt.loglog(ksub,pksub.real*ratio,label='single subtracting host')

file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/%s_singlewohostp_%s_%s.txt'% (fit_lab,name,numb), 'w')
for i in pksub.real:
    file1.write('%s\n' % i)
file1.close()
file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/%s_singlewohostk_%s_%s.txt' % (fit_lab,name,numb), 'w')
for i in ksub:
    file1.write('%s\n' % i)
file1.close()"""

#################################################################
#average subtracting the host
#################################################################

allspectra = []
allks = []
for i in allprojs:
    sub = i - host 
    sub_mesh = ArrayMesh(sub,(bs,bs),Nmesh=sub.shape[0])
    Pk2 = ProjectedFFTPower(sub_mesh,axes=(0,1),BoxSize=bs,Nmesh=sub.shape[0])
    allspectra.append(Pk2.power['power'])
    allks.append(Pk2.power['k'])

ps2 = np.mean(allspectra,axis=0)
ks2 = np.mean(allks,axis=0)

file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/avgwohostp_%s_%s.txt'% (name,numb), 'w')
for i in ps2.real:
    file1.write('%s\n' % i)
file1.close()
file1 = open('/n/home04/adiazrivero/conv_powerspec/particle_level/powerspectra/avgwohostk_%s_%s.txt'% (name,numb), 'w')
for i in ks2:
    file1.write('%s\n' % i)
file1.close()


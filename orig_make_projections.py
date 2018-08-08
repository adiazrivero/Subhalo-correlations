from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from functions import *
import cPickle as pickle
from tburk_convergence import *
from tnfw_convergence import *

start_time = time.time()

#executing variables declared externally

list=sys.argv[1:]
for i in list:
    exec(i)

if name == 0:
    name = 'CDM'
elif name == 1:
    name = 'ETHOS_1'
elif name == 2:
    name = 'ETHOS_4'

print name

################################################

if numb == 78 or numb == 95:
    numb = '0' + str(numb)
else:
    numb = str(numb)

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

################################################

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

################################################

list1 = []
for i in range(16):
    list1.append('/n/hernquistfs3/jzavala/ETHOS/%s/groups_%s/fof_subhalo_tab_%s.%s.hdf5' % (name,numb,numb,i))    

#obtaining subhalo catalogs

mpart = 1.90382332e-6 / 0.6909
m_min = 50 * mpart

m_min = m_min * 1e10
#m_high = 1e15
m_high = mhigh * 1e10

print m_min,m_high

mass,pos,rh,rvmax = simul_data(list1,mlow=m_min,mhigh=m_high,Rnge=600) #mhigh is in units of 10^10 M_sun

print "total number of subhalos within L_box < 300 kpc: %s " % np.shape(mass)
print "max mass = %e" % np.amax(mass)
print "min mass = %e" %  np.amin(mass)

"""mass,pos,rh,rvmax = read_catalogue(list1, mlow=m_min, mhigh=m_high,L_box=600)

print "total number of subhalos within L_box < 300 kpc: %s " % np.shape(mass)
print "max mass = %e" % np.amax(mass)
print "min mass = %e" %  np.amin(mass)

sys.exit()"""

dicti = pickle.load(open('/n/home04/adiazrivero/conv_powerspec/catalog_level/host_features/host_%s_%s.txt' % (name,numb),'rb'))

Mhost = dicti['host_mass']

arr = np.load('/n/home04/adiazrivero/conv_powerspec/catalog_level/rein_sigmac_all/rein_sigmac_%s_%s.npy' % (name,numb))
zs = arr[0]
zl = arr[1]
rein = arr[3]
sigmac = arr[4]

if zl != z:
    print "The wrong simulation is being used as a lens!"
    sys.exit()

print "r_ein = %s kpc " % rein
print "sigma_crit = %e M_sun/kpc^2 " % sigmac

shft = 0
rnge = 600

posmass = orig_projection(pos,mass,rh,rvmax,rnge=rnge,shift=shft,num_proj=1)

posmass2 = np.asarray(posmass)
print "proj num = %s " % len(posmass2)

dir1 = '/n/home04/adiazrivero/conv_powerspec/catalog_level/%s/%s/%s/origproj_%s_%s_%s' % (name,numb,mlab,name,numb,mlab)
np.save(dir1, posmass2)

print 'original projection obtained!'
print 'saved in %s' % dir1


arr = np.load('/n/home04/adiazrivero/conv_powerspec/catalog_level/rein_sigmac_all/rein_sigmac_%s_%s.npy' % (name,numb))

zs = arr[0]
zl = arr[1]
rein = arr[3]
sigmac = arr[4]

if zl != z:
    print "The wrong simulation is being used as a lens!"
    sys.exit()

pos_2d = []
masses = []
rmaxs = []
rhs = []
for i in posmass:
    r,m,rh,rmax,_ = zip(*i)
    pos_2d.append(r)
    masses.append(m)
    rmaxs.append(rmax)
    rhs.append(rh)

masses = np.array([np.array(i) for i in masses])
rhs = np.array([np.array(i) for i in rhs])
rmaxs = np.array([np.array(i) for i in rmaxs])

taus = []
mscale = []

if name == 'CDM':

    tau_vs_xhalf_arr = np.loadtxt('tau_vs_xhalf.txt')
    tau_vs_xhalf = interp1d(tau_vs_xhalf_arr[:,0],tau_vs_xhalf_arr[:,1])

    r_s = rmaxs/2.1626


    for i,j,k in zip(rhs,r_s,masses):

        tau = tau_vs_xhalf(i/j)
        m_NFW = k * ((tau**2 + 1)**2/(tau**2*((tau**2 - 1)*np.log(tau) + np.pi*tau - (tau**2 + 1))))

        taus.append(np.array(tau))
        mscale.append(np.array(m_NFW))

elif name == 'ETHOS_4':

    tau_vs_xhalf_arr = np.loadtxt('tau_vs_xhalf_burk.txt')
    tau_vs_xhalf = interp1d(tau_vs_xhalf_arr[:,0],tau_vs_xhalf_arr[:,1])

    r_b = rmaxs/3.3446
    p = 0.66651

    for i,j,k in zip(rhs,r_b,masses):

        select1 = (i/j) < min(tau_vs_xhalf_arr[:,0])
        select2 = (i/j) > max(tau_vs_xhalf_arr[:,0])
        #print len(i/j),len((i/j)[select1 ^ select2]),len((i/j)[~(select1 ^ select2)]) 

        tau = np.zeros( (i/j).shape )
        tau[~(select1 ^ select2)] = tau_vs_xhalf( (i/j)[~(select1 ^ select2)] )

        ex = extrap1d(tau_vs_xhalf)
        tau[(select1 ^ select2)] = ex( (i/j)[(select1 ^ select2)] )

        #print tau

        n = (tau**2 * (np.pi * (p-tau)**2 + 4*tau**2 * np.log(p/tau)))/(4 * (p**4 - tau**4))
        m_burk = k / n

        taus.append(np.array(tau))
        mscale.append(np.array(m_burk))

    r_s = r_b / p

r_t = np.array(taus) * r_s

pix_num = 501
pix_size = rnge/pix_num
x = np.linspace(-rnge/2, rnge/2, pix_num)
y = np.linspace(-rnge/2, rnge/2, pix_num)

count = 0
conv_list = []

"""r_s = 63/2.1626
#r_s = 4
t = 244/r_s
#t = 18.2784438688
m_nfw = mnfw(1.6e12,t)
print 'mnfw= %e' % m_nfw
result = conv_tnfw(x[:,None]-0.1,y[None,:]-0.1,t,m_nfw,r_s).real
print 'max kappa = %.3f' % np.amax(result)
sys.exit()"""

for g,q,u,z in zip(pos_2d,taus,mscale,r_s):
    count += 1
    print count
    beginning = time.time()

    if name == 'CDM':
        print 'using the truncated NFW profile'
        result = [conv_tnfw(x[:,None]-i[0],y[None,:]-i[1],k,l,m,sigmac).real for i,k,l,m in zip(g,q,u,z)]

    elif name == 'ETHOS_4':
        print 'using the truncated Burkert profile'
        #def conv_tburk(x,y,rs,mb,pval,tauval,sigmac):        
        result = [conv_tburk(x[:,None]-i[0],y[None,:]-i[1],m,l,p,k,sigmac).real for i,k,l,m in zip(g,q,u,z)]

    res = sum(result)
    print ' ----------------------- '
    conv_list.append(res)
    print("%s seconds" % (time.time() - beginning))

    conv_list2 = np.asarray(conv_list)

    dir2 = '/n/home04/adiazrivero/conv_powerspec/catalog_level/%s/%s/%s/origconv%s_%s_%s_%s_%s_%s_2' % (name,numb,mlab,name,numb,pix_num,count,rnge,mlab)
    #dir2 = '/n/home04/adiazrivero/conv_powerspec/catalog_level/%s/%s/%s/origconv%s_%s_%s_%s_%s_%s' % (name,numb,mlab,name,numb,pix_num,count,rnge,mlab)
    np.save(dir2,conv_list2)
    print 'saved in %s' % dir2


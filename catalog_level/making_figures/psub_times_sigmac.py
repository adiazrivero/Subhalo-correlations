from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys

plt.rc('text', usetex=True)

bs = 600
h = 0.6909

list=sys.argv[1:]
for i in list:
    exec(i)

if name == 0:
    name = 'CDM'
    fit_lab = 'tnfw'
elif name == 1:
    name = 'ETHOS_1'
elif name == 2:
    name = 'ETHOS_4'
    fit_lab = 'tburk'

print name

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

mhigh = mhigh * 1e10

######################################################

mlab2 = 'm6'
rnge = 100

if name == 'ETHOS_4':

    cls = ['green','cyan','orange']

elif name == 'CDM':

    cls = ['green','blue','orange']

for numb,clr in zip(['127','095','078'],cls):

    arr = np.load('/Users/anadiazrivero/Dropbox/Research/Gravitational_lensing/Simulations_paper2/subhalo_catalogs/%s/%s/%s/ind_curves_%s_%s_501_%s_%s.npy' % (name,numb,mlab,name,numb,rnge,mlab2))

    sigmac = np.load('/Users/anadiazrivero/Dropbox/Research/Gravitational_lensing/Simulations_paper2/rein_sigmac_all/rein_sigmac_%s_%s.npy' % (name,numb))[-1]

    p1 = np.percentile(arr,2.5,axis=0)
    p2 = np.percentile(arr,16,axis=0)
    p3 = np.percentile(arr,50,axis=0)
    p4 = np.percentile(arr,84,axis=0)
    p5 = np.percentile(arr,97.5,axis=0)

    ks = np.loadtxt('/Users/anadiazrivero/Dropbox/Research/Gravitational_lensing/Simulations_paper2/subhalo_catalogs/%s/%s/%s/k%s_%s_501_%s_%s.txt'% (name,numb,mlab,name,numb,rnge,mlab2))

    if numb == '127':
        label = '$z=0^*$'
    elif numb == '095':
        label = '$z=0.5$'
    elif numb == '078':
        label = '$z=1$'

    plt.loglog(ks,p3,color=clr,label='%s' % label)
    plt.fill_between(ks,p2,p4,alpha=0.3,color=clr)
    plt.fill_between(ks,p1,p5,alpha=0.1,color=clr)

    plt.xlim(np.amin(ks),10)
    plt.tick_params(axis='both', which='major', labelsize=35)
    plt.xlabel('$k$ [kpc$^{-1}$]',fontsize=35)
    #plt.ylabel('$P_{\\rm sub}(k) \\Sigma_{\\rm crit}$ [kpc$^2$]',fontsize=35)
    plt.ylabel('$P_{\\rm sub}(k) \\Sigma_{\\rm crit}$ [M$_{\\odot}$]',fontsize=35)

if name == 'ETHOS_4':

    plt.text(7*0.01,1.5*1e4,'ETHOS4',fontsize=40)

elif name == 'CDM':

    plt.text(7*0.01,1.5*1e4,'CDM',fontsize=40)

plt.legend(prop={'size': 30})

plt.ylim(1e0,4*1e4)
plt.gcf().subplots_adjust(bottom=0.17)
plt.show()

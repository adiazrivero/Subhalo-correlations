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

######################################################

mlab2 = 'm6'
rnge=100

for mlab,clr in zip(['m11','m9','m8'],['red','purple','cyan']):

    arr = np.load('/Users/anadiazrivero/Dropbox/Research/Gravitational_lensing/Simulations_paper2/subhalo_catalogs/%s/%s/%s/ind_curves_%s_%s_501_%s_%s.npy' % (name,numb,mlab,name,numb,rnge,mlab2))

    p1 = np.percentile(arr,2.5,axis=0)
    p2 = np.percentile(arr,16,axis=0)
    p3 = np.percentile(arr,50,axis=0)
    p4 = np.percentile(arr,84,axis=0)
    p5 = np.percentile(arr,97.5,axis=0)

    ks = np.loadtxt('/Users/anadiazrivero/Dropbox/Research/Gravitational_lensing/Simulations_paper2/subhalo_catalogs/%s/%s/%s/k%s_%s_501_%s_%s.txt'% (name,numb,mlab,name,numb,rnge,mlab2))

    if mlab == 'm11':
        label = 'All subhalos'
    elif mlab == 'm9':
        label = '$m_{\\rm high} = 10^9$ M$_{\\odot}$'
    elif mlab == 'm8':
        label = '$m_{\\rm high} = 10^8$ M$_{\\odot}$'

    plt.loglog(ks,p3,color=clr,label='%s' % label)
    plt.fill_between(ks,p2,p4,alpha=0.3,color=clr)
    plt.fill_between(ks,p1,p5,alpha=0.1,color=clr)

    plt.xlim(np.amin(ks),10)
    plt.ylim(1e-9,1e-2)
    plt.tick_params(axis='both', which='major', labelsize=35)
    plt.xlabel('$k$ [kpc$^{-1}$]',fontsize=35)
    plt.ylabel('$P_{\\rm sub}(k)$ [kpc$^2$]',fontsize=35)

    if name == 'ETHOS_4':
        plt.text(7*0.01,2.5*1e-3,'ETHOS4',fontsize=40)
    elif name == 'CDM':
        plt.text(7*0.01,2.5*1e-3,'CDM',fontsize=40)

    plt.legend(prop={'size': 30})

plt.ylim(1e-9,1e-5)
plt.gcf().subplots_adjust(bottom=0.17)
plt.show()

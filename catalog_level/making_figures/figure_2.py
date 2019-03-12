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
numb='095'
rnge=200

fig,ax = plt.subplots(1)

for name,clr,label in zip(['CDM','ETHOS_4'],['blue','cyan'],['CDM','ETHOS4']):

    arr = np.load('/Users/anadiazrivero/Dropbox/convergence_ps_paper2/clean_fof_masscut/%s/%s/%s/ind_curves_%s_%s_501_%s_%s.npy' % (name,numb,mlab,name,numb,rnge,mlab2))
    ks = np.loadtxt('/Users/anadiazrivero/Dropbox/convergence_ps_paper2/clean_fof_masscut/%s/%s/%s/k%s_%s_501_%s_%s.txt'% (name,numb,mlab,name,numb,rnge,mlab2))

    p50 = np.percentile(arr,50,axis=0)
    p68 = np.percentile(arr,68,axis=0)
    p32 = np.percentile(arr,32,axis=0)
    p95 = np.percentile(arr,95,axis=0)
    p5 = np.percentile(arr,5,axis=0)

    if name == 'CDM':
        #plus = (1.775e-6 + 2.626e-6) * np.ones(len(ks))
        plus = (1.675e-6 + 9.239e-7) * np.ones(len(ks))
        minus = (1.675e-6 - 9.239e-7) * np.ones(len(ks))

        ax.fill_between(ks,minus,plus,color='gray',alpha=0.3)
        ax.text(0.75,1e-6,'$\\bar{\\kappa}_{\\rm sub} m_{\\rm eff} / \\Sigma_{\\rm crit}|_{\\rm CDM}$',fontsize=40)

        k_s = 1 / 0.0493629205658
        k_t = 1 / 12.2838205247

        print k_s,k_t

        ax.axvline(k_s,color='k',linestyle='dotted')
        ax.axvline(k_t,color='k',linestyle='dashed',label='$k_{\\rm trunc, CDM}$')

        nbar2sh = np.loadtxt('/Users/anadiazrivero/Dropbox/convergence_ps_paper2/nbar_2D_2sh_%s_%s_501_%s_%s.txt' % (name,numb,rnge,mlab2))
        ax.loglog(ks,nbar2sh,color='red',linestyle='dashed',label='$P_{\\rm 2sh,CDM}(k)$')

    ax.loglog(ks,p50,color=clr,label='%s' % label)
    ax.fill_between(ks,p32,p68,alpha=0.3,color=clr) #,label='68th perc')
    ax.fill_between(ks,p5,p95,alpha=0.1,color=clr) #,label='95th perc')

    ax.set_xlim(np.amin(ks),6)
    ax.set_ylim(1e-9,1e-5)
    ax.tick_params(axis='both', which='major', labelsize=35)
    ax.set_xlabel('$k$ [kpc$^{-1}$]',fontsize=35)
    ax.set_ylabel('$P_{\\rm sub}(k)$ [kpc$^2$]',fontsize=35)

    if name == 'ETHOS_4':
        #plt.text(7*0.01,4.5*1e-6,'ETHOS4',fontsize=40)
        plt.text(7*0.01,2.5*1e-3,'ETHOS4',fontsize=40)
    elif name == 'CDM':
        #plt.text(7*0.01,4.5*1e-6,'CDM',fontsize=40)
        plt.text(7*0.01,2.5*1e-3,'CDM',fontsize=40)

ax.legend(prop={'size': 30},loc=3)
#plt.savefig('allcurves_%s.png' % mlab)

plt.gcf().subplots_adjust(bottom=0.17)
#plt.tight_layout()
plt.show()


sys.exit()

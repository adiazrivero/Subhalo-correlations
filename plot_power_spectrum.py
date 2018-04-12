from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import cmath
from R_ein import rein_sigmac
from functions import simul_data,rotation,projections,twoD_ps,angular_average,variance

np.set_printoptions(suppress=True,linewidth=np.nan,threshold=np.nan)
plt.rc('text', usetex=True)

plot_errors = True
save_files = False

############################################
#loading properties of the convergence maps
############################################

feats = np.load('convfeat_101_100.npy').item()
num_proj = feats.get('num_proj')
rnge = feats.get('rnge')
pix_num = feats.get('pixels')
shft = feats.get('shft')
m_high_cut = feats.get('m_high cut')
pix_size = rnge/pix_num

projfeats = np.load('projfeat_%s_%s.npy' % (pix_num,rnge)).item()
rein = projfeats.get('r_ein')
sigmac = projfeats.get('sigmac')
meff = projfeats.get('m_eff')
mavg = projfeats.get('m_avg')
avgkappa = projfeats.get('k_avg')
avgkappa_error = projfeats.get('k_avg_error')
max_rt = projfeats.get('avg_max_rt')
min_rt = projfeats.get('avg_min_rs')

############################################
#loading the convergence maps
############################################

conv_list1 = np.load('conv_%s_10_%s.npy' % (pix_num,rnge))
conv_list2 = np.load('conv_%s_20_%s.npy' % (pix_num,rnge))
conv_list3 = np.load('conv_%s_30_%s.npy' % (pix_num,rnge))

cl1 = np.ndarray.tolist(conv_list1)
cl2 = np.ndarray.tolist(conv_list2)
cl3 = np.ndarray.tolist(conv_list3)

conv_list = []
for i,j,k in zip(cl1,cl2,cl3):
    conv_list.append(i)
    conv_list.append(j)
    conv_list.append(k)

avg_conv = np.mean(conv_list,axis=0)

####################################################################
# Doing the 2D FT and angular-averaging to obtain the PS
####################################################################

individual_ps,tot_ps,kx,ky = twoD_ps(data=conv_list,pix_size=pix_size,rnge=rnge,shift=shft,show_ps=False)

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
ps1d,K = angular_average(tot_ps,kx,ky,mask=mask,rnge=rnge,pix_num=pix_num,dr=pix_size_k,remove_first=True)

amp = avgkappa * meff / sigmac
print 'PS amplitude = %s ' % amp

plt.title('Convergence field on the lens plane (%s proj/$\kappa_{avg}$=%.5f/$m_{high}$ cut = %s)' % (len(conv_list),avgkappa,m_high_cut))
ax = plt.subplot(111)
ax.axhline(amp,0.01,100,c='c')
ax.set_xscale("log")
ax.set_yscale("log")

if plot_errors == True:
    var = variance(individual_ps,ps1d,len(individual_ps),kx,ky,mask=mask,rnge=rnge,pix_num=pix_num,pix_size=pix_size_k)
    std = [np.sqrt(i) for i in var]
    min_err = [i-j for i,j in zip(ps1d,std)]
    plus_err = [i+j for i,j in zip(ps1d,std)]

    ax.plot(K,ps1d)
    ax.fill_between(K,min_err,plus_err,alpha=0.4,facecolor='red')

elif std == False:
    print 'b'
    ax.plot(K,ps1d)

ax.set_xlabel('$k$ [kpc]$^{-1}$')
ax.set_ylabel('$P(k)$ [kpc]$^2$')
#plt.legend()
plt.show()

if save_files == True:

    file1 = open('ps_%s_%s.txt'% (pix_num,rnge), 'w')
    for i in ps1d:
        file1.write('%s\n' % i)
    file1.close()

    file2 = open('k_%s_%s.txt'% (pix_num,rnge), 'w')
    for i in K:
        file2.write('%s\n' % i)
    file2.close()

    file3 = open('min_%s_%s.txt'% (pix_num,rnge), 'w')
    for i in min_err:
        file3.write('%s\n' % i)
    file3.close()

    file4 = open('plus_%s_%s.txt'% (pix_num,rnge), 'w')
    for i in plus_err:
        file4.write('%s\n' % i)
    file4.close()

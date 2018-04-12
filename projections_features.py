from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pylab as py
import sys
import time
#from rein import rein_sigmac
from functions import simul_data,rotation,projections,angular_average,twoD_ps

plt.rc('text', usetex=True)

############### loading data #####################

posmass = np.load('projections.npy')
m_high_cut = False

feats = np.load('convfeat_101_100.npy').item()

num_proj = feats.get('num_proj')
rnge = feats.get('rnge')
pix_num = feats.get('pixels')
shft = feats.get('shft')
m_high_cut = feats.get('m_high cut')
max_rt = feats.get('avg_max_rt')

pix_size = rnge/pix_num

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

##################################################

rein = 10.300324225
sigmac = 3120194565.15

pos_2d = []
masses = []
for i in posmass:
    r,m,_,_ = zip(*i)
    pos_2d.append(r)
    masses.append(m)

num_sr = []
for i in pos_2d:
    subh_rein = [j for j in i if  np.sqrt(j[0]**2+j[1]**2) <= rein]
    num_sr.append(len(subh_rein))

print 'average num of subh in SL region = %s +- %s ' % (np.mean(num_sr),np.std(num_sr))

num = [len(i) for i in masses]
mean_mass = [np.mean(i) for i in masses]

print 'average num of subh r < 100 kpc = %s ' % np.mean(num),np.std(num)
print 'average m_sub = %s' % np.mean(mean_mass),np.std(mean_mass)

m2av = [j**2 for j in i for i in masses]

m2avg = np.asarray([np.mean(i) for i in m2av])
m2 = np.mean(m2avg)
meff = np.exp(np.log(m2)-np.log(np.mean(mean_mass)))

print "m_eff = %s " % meff

def rs(msub):
    return 0.1*(msub/10**6)**(1/3)

def mnfw(m,tau):
    n = np.asarray([(i**2/(i**2 + 1)**2) * ((i**2 - 1) * np.log(i) + i * np.pi - (i**2 + 1)) for i in tau])
    mnfw = m/n
    return mnfw

masses = np.asarray([np.asarray(i) for i in masses])
r_s = rs(masses)
avg_min_rs = np.mean([min(i) for i in r_s])
rt = 5*r_s
max_rt = [max(i) for i in rt]
avg_max_rt = np.mean(max_rt)

print 'average r_s,min = %s +- %s ' % (avg_min_rs,np.std([min(i) for i in r_s]))
print 'average r_t,max = %s +- %s ' % (avg_max_rt,np.std(max_rt))

x = np.linspace(-rnge/2, rnge/2, pix_num)
y = np.linspace(-rnge/2, rnge/2, pix_num)

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

kappar,r = angular_average(avg_conv,x,y,mask=mask,rnge=rnge,pix_num=pix_num,dr=pix_size,remove_first=False)
print "Average convergence in the SL region: %s +- %s " % (np.mean(kappar),np.std(kappar))

#plotting the average convergence map
fig,(ax1) = plt.subplots(1)
k = ax1.imshow(np.log10(avg_conv),extent=[-rnge/2,rnge/2,-rnge/2,rnge/2],origin='lower')
circ = Circle((0,0),rein,facecolor='none',fill=False,linestyle='dashed')
ax1.add_patch(circ)
fig.colorbar(k)
ax1.set_ylabel('$kpc$',fontsize=15)
ax1.set_xlabel('$kpc$',fontsize=15)
ax1.set_title('Convergence field on the lens plane (%s proj/$\\bar{\kappa}_{sub}$=%.6f/$m_{high}$ cut=%s)' % (len(posmass),np.mean(kappar),m_high_cut),fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.show()

projfeat = {'r_ein':rein,'sigmac':sigmac,'m_high cut':m_high_cut,'m_avg':np.mean(mean_mass),'m_eff': meff,'k_avg':np.mean(kappar),'k_avg_error':np.std(kappar),"avg_max_rt":avg_max_rt,"avg_min_rs":avg_min_rs}

np.save('projfeat_%s_%s' % (pix_num,rnge), projfeat)

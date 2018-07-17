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

if mhigh == 0:
   mlab = 'm0'
if mhigh == 1e-2:
   mlab = 'm8'
if mhigh == 1e-3:
   mlab = 'm7'

print numb,name,mlab

posmass = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/projections_%s_%s_%s.npy' % (name,numb,mlab,name,numb,mlab))
pix_num0 = 501
rnge0 = 100


if name == 'CDM':
    #feats = np.load('convfeat%s_%s_501_100.npy' % (name,numb)).item()
    feats = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convfeat%s_%s_%s_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num0,rnge0,mlab)).item()
elif name == 'ETHOS_1' or name == 'ETHOS_4':
    #feats = np.load('convburkfeat%s_%s_501_100.npy' % (name,numb)).item()
    feats = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburkfeat%s_%s_%s_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num0,rnge0,mlab)).item()


num_proj = feats.get('num_proj')
rnge = feats.get('rnge')
pix_num = feats.get('pixels')
shft = feats.get('shft')
m_high_cut = feats.get('m_high cut')
max_rt = feats.get('avg_max_rt')

if pix_num0 != pix_num or rnge0 != rnge:
    print "Using wrong projections!"
    sys.exit()

pix_size = rnge/pix_num

if name == 'CDM':
   
    conv_list1 = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/conv%s_%s_%s_10_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num,rnge,mlab))    
    conv_list2 = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/conv%s_%s_%s_20_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num,rnge,mlab))    
    conv_list3 = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/conv%s_%s_%s_30_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num,rnge,mlab))    

elif name == 'ETHOS_1' or name == 'ETHOS_4':

    conv_list1 = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburk%s_%s_%s_10_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num,rnge,mlab))    
    conv_list2 = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburk%s_%s_%s_20_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num,rnge,mlab))    
    conv_list3 = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburk%s_%s_%s_30_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num,rnge,mlab))    

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

arr = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/rein_sigmac_all/rein_sigmac_%s_%s.npy' % (name,numb))

zs = arr[0]
zl = arr[1]
rein = arr[3]
sigmac = arr[4]

if zl != z:
    print "The wrong simulation is being used as a lens!"
    sys.exit()

pos_2d = []
masses = []
for i in posmass:
    r,m,_,_,_ = zip(*i)
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
m2err = np.std(m2avg)
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

"""print 'new k_sub = %s' % np.mean(avg_conv)
print np.std(avg_conv)
sys.exit()"""

#plotting the average convergence map
"""fig,(ax1) = plt.subplots(1)
k = ax1.imshow(np.log10(avg_conv),extent=[-rnge/2,rnge/2,-rnge/2,rnge/2],origin='lower')
circ = Circle((0,0),rein,facecolor='none',fill=False,linestyle='dashed')
ax1.add_patch(circ)
fig.colorbar(k)
ax1.set_ylabel('$kpc$',fontsize=15)
ax1.set_xlabel('$kpc$',fontsize=15)
ax1.set_title('Convergence field on the lens plane (%s proj/$\\bar{\kappa}_{sub}$=%.6f/$m_{high}$ cut=%s)' % (len(posmass),np.mean(kappar),m_high_cut),fontsize=15)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.show()"""

projfeat = {'m^2': m2, 'sigma(m^2)': m2err,'r_ein':rein,'sigmac':sigmac,'m_high cut':m_high_cut,'m_high':mhigh,'m_avg':np.mean(mean_mass),'sigma(m_avg)':np.std(mean_mass),'m_eff': meff,'k_avg':np.mean(avg_conv),'k_avg_error':np.std(avg_conv),"avg_max_rt":avg_max_rt,'sigma(r_t,max)':np.std(max_rt),"avg_min_rs":avg_min_rs,'sigma(r_s,min)':np.std([min(i) for i in r_s]),'N_rein': np.mean(num_sr), 'sigma(N_rein)':np.std(num_sr), 'N_100': np.mean(num),'sigma(N_100)':np.std(num)}
#projfeat = {'m^2': m2, 'sigma(m^2)': m2err,'r_ein':rein,'sigmac':sigmac,'m_high cut':m_high_cut,'m_high':mhigh,'m_avg':np.mean(mean_mass),'sigma(m_avg)':np.std(mean_mass),'m_eff': meff,'k_avg':np.mean(kappar),'k_avg_error':np.std(kappar),"avg_max_rt":avg_max_rt,'sigma(r_t,max)':np.std(max_rt),"avg_min_rs":avg_min_rs,'sigma(r_s,min)':np.std([min(i) for i in r_s]),'N_rein': np.mean(num_sr), 'sigma(N_rein)':np.std(num_sr), 'N_100': np.mean(num),'sigma(N_100)':np.std(num)}

if name == 'CDM':
    #np.save('projfeat%s_%s_%s_%s' % (name,numb,pix_num,rnge), projfeat)
    np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/projfeat%s_%s_%s_%s_%s' % (name,numb,mlab,name,numb,pix_num,rnge,mlab),projfeat)

elif name == 'ETHOS_1' or name == 'ETHOS_4':
    np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/projfeatburk%s_%s_%s_%s_%s' % (name,numb,mlab,name,numb,pix_num,rnge,mlab),projfeat)
    #np.save('projfeatburk%s_%s_%s_%s' % (name,numb,pix_num,rnge), projfeat)






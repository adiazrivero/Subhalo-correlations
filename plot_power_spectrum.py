from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import cmath
from functions import *

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

pix_num0 = 501
rnge0 = 100

############################################
#loading properties of the convergence maps
############################################

if name == 'CDM':
    
    feats = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convfeat%s_%s_%s_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num0,rnge0,mlab)).item()
    projfeats = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/projfeat%s_%s_%s_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num0,rnge0,mlab)).item()

elif name == 'ETHOS_1' or name == 'ETHOS_4':

    feats = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburkfeat%s_%s_%s_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num0,rnge0,mlab)).item()
    projfeats = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/projfeatburk%s_%s_%s_%s_%s.npy' % (name,numb,mlab,name,numb,pix_num0,rnge0,mlab)).item()


num_proj = feats.get('num_proj')
rnge = feats.get('rnge')
pix_num = feats.get('pixels')
shft = feats.get('shft')
m_high_cut = feats.get('m_high cut')
pix_size = rnge/pix_num

if pix_num0 != pix_num or rnge0 != rnge:
    print "Using wrong projections!"
    sys.exit()

rein = projfeats.get('r_ein')
sigmac = projfeats.get('sigmac')
meff = projfeats.get('m_eff')
mavg = projfeats.get('m_avg')
avgkappa = projfeats.get('k_avg')
avgkappa_error = projfeats.get('k_avg_error')
max_rt = projfeats.get('avg_max_rt')
min_rt = projfeats.get('avg_min_rs')


amp = avgkappa * meff / sigmac
print 'PS amplitude = %s ' % amp

############################################
#loading the convergence maps
############################################

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

#obtain the 1d power spectrum (for any n) and errors

pix_size_k = np.abs(kx[0]-kx[1])

ns = [0,1,2]
power_spectra,K = multipoles(tot_ps,kx,ky,mask=mask,pix_num=pix_num,dr=pix_size_k,ns=ns)

K = K[1:]
    
file0 = open('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/k%s_%s_%s_%s.txt'% (name,numb,mlab,name,numb,pix_num,rnge), 'w')
for i in K:
    file0.write('%s\n' % i)
file0.close()

for key in power_spectra.keys():
 
    var = variance(individual_ps,power_spectra[key][1:],len(individual_ps),kx,ky,mask=mask,rnge=rnge,pix_num=pix_num,pix_size=pix_size_k,n=key)
    std = [np.sqrt(i) for i in var]
    min_err = [i-j for i,j in zip(power_spectra[key][1:],std)]
    plus_err = [i+j for i,j in zip(power_spectra[key][1:],std)]

    if save_files == True:

        if key == '0':
	    lab = 'monopole'
	elif key == '1':
	    lab = 'dipole'
	elif key == '2':
	    lab = 'quadrupole'

        file1 = open('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/%s%s_%s_%s_%s.txt'% (name,numb,mlab,lab,name,numb,pix_num,rnge), 'w')
        for i in power_spectra[key][1:]:
            file1.write('%s\n' % i)
        file1.close()
    
        file2 = open('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/min_%s%s_%s_%s_%s.txt'% (name,numb,mlab,lab,name,numb,pix_num,rnge), 'w')
        for i in min_err:
            file2.write('%s\n' % i)
        file2.close()

        file3 = open('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/plus_%s%s_%s_%s_%s.txt'% (name,numb,mlab,lab,name,numb,pix_num,rnge), 'w')
        for i in plus_err:
            file3.write('%s\n' % i)
        file3.close()
    


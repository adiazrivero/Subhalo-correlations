from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import h5py
import sys
import itertools
import scipy
import time 
import json
import cPickle as pickle

def host_data(data):
    masses = []
    positions1 = []
    rh = []
    rvmax = []
    for i in data:
        file = h5py.File(i,'r')
        h = file['Header'].attrs['HubbleParam']
        
	if len(file['Subhalo'].keys()) != 0:
            masses.append(file['Subhalo']['SubhaloMass'].value) #in 10^10 M_sun/h
            positions1.append(file['Subhalo']['SubhaloPos'].value) # in kpc/h, absolute box coordinates
            rh.append(file['Subhalo']['SubhaloHalfmassRad'].value) # in kpc/h
            rvmax.append(file['Subhalo']['SubhaloVmaxRad'].value) # in kpc/h
    
    print 'redshift = %s' % file['Header'].attrs['Redshift']
    print 'h = %s' % h

    masses = [i/h for j in masses for i in j]
    positions1 = [np.asarray(i)/h for j in positions1 for i in j]
    rh = [i/h for j in rh for i in j]
    rvmax = [i/h for j in rvmax for i in j]
    
    parent_mass = masses[0] * 10**(10)
    parent_pos = positions1[0]
    parent_rvmax = rvmax[0]
    parent_rh = rh[0]
    
    return {'host_mass':parent_mass,'host_center':parent_pos,'host_rmax':parent_rvmax,'host_rh':parent_rh}

start_time = time.time()

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
    print "z = 0"
elif numb == '095':
    print 'z = 0.5'
elif numb == '078':
    print 'z = 1'
else:
    print "wrong redshift specified"
    sys.exit()

list1 = ['/n/hernquistfs3/jzavala/ETHOS/%s/groups_%s/fof_subhalo_tab_%s.0.hdf5' % (name,numb,numb)]
host_dict = host_data(list1)

print "host mass = %e" % host_dict['host_mass']
#print "host R_max = %s" % host_dict['host_rmax']
#print "host R_1/2 = %s" % host_dict['host_rh']
print 'host center = %s' % host_dict['host_center']

#with open('host_%s_%s.txt' % (name,numb), 'w') as file:
     #file.write(json.dumps(host_dict))

with open('host_%s_%s.txt' % (name,numb), 'w') as file:
     file.write(pickle.dumps(host_dict))


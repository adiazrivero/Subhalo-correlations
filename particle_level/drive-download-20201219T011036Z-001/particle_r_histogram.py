import numpy as np
import matplotlib as mpl
#import pylab as pl
import h5py
import sys
import time
#from mpl_toolkits.mplot3d import Axes3D

#######################################################################
#importing data
#######################################################################

start_time = time.time()


num = '127'
list = []
for i in range(16):
    list.append('/n/hernquistfs3/jzavala/ETHOS/CDM/snapdir_%s/snap_%s.%s.hdf5' % (num,num,i))

coords = []
for i in list:
    file = h5py.File(i,'r')
    for i in file.keys():
	print i
    
    for i in file['Header'].attrs:
        print i
    
    pm = file['Header'].attrs['MassTable'][1] * 1e10
    print 'particle mass = %s' % pm
    
    print 'redshift = %s' % file['Header'].attrs['Redshift']
    print 'box size = %s' % file['Header'].attrs['BoxSize']
    print 'hubble parameter = %s' % file['Header'].attrs['HubbleParam']
    
    coords2 = file['PartType2']['Coordinates'].value
    masses = file['PartType2']['Masses'].value
    coords.append(coords2)


coord = [i[0:2] for j in coords for i in j]
x,y = zip(*coord)

r = [np.sqrt(i**2+j**2) for i,j in zip(x,y)]
plt.hist(r,bins=50)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

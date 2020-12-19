from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
import h5py
from nbodykit.source.catalog import HDFCatalog,CSVCatalog
from nbodykit import transform
from nbodykit.lab import ProjectedFFTPower, FFTPower,BigFileMesh
from nbodykit.source.mesh.array import ArrayMesh
from rotation import rotation
import dask.array as da
import cPickle as pickle

start_time = time.time()

inputs = sys.argv[1:]
for i in inputs:
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

lim = 300
bs = 2 * lim
nmesh = 1400
h = 0.6909

dicti = pickle.load(open('/n/home04/adiazrivero/conv_powerspec/catalog_level/host_%s_%s.txt' % (name,numb),'rb'))
shiftx = dicti['host_center'][0]
shifty = dicti['host_center'][1]
shiftz = dicti['host_center'][2]

f = CSVCatalog('/n/dvorkin_lab/anadr/%s_catalog/%s/%s_%s_parttype1_*.txt' % (name,numb,name,numb), ['x','y','z'])

f['x'] = f['x'] / h - shiftx
f['y'] = f['y'] / h - shifty
f['z'] = f['z'] / h - shiftz


select1 = f['x'] < lim
select2 = f['x'] > -lim
select3 = f['y'] < lim
select4 = f['y'] > -lim
select5 = f['z'] < lim
select6 = f['z'] > -lim

subcat = f[select1 & select2 & select3 & select4 & select5 & select6]
subcat['Position'] =  transform.StackColumns(subcat['x'], subcat['y'], subcat['z'])

mesh = subcat.to_mesh(Nmesh=nmesh,BoxSize=bs)
mesh.save('/n/dvorkin_lab/anadr/%s_%s_%s_%s_hbugfix.bigfile' % (name,numb,nmesh,bs))

N_part = len(subcat)
print 'number of particles: %s' % N_part

box_vol = bs**3
nbar = N_part / box_vol

print nbar

dicti = {'nbar': nbar, 'N_part': N_part}

with open('/n/home04/adiazrivero/conv_powerspec/particle_level/nbar/hbugfix/nbar_%s_%s_%s_%s.txt' % (name,numb,nmesh,bs), 'w') as file:
     file.write(pickle.dumps(dicti))

print("%s seconds" % (time.time() - start_time))


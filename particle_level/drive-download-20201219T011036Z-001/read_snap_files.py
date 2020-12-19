import numpy as np
import h5py 
import matplotlib.pyplot as plt 
import sys 

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

for k in range(16):
    file = h5py.File('/n/hernquistfs3/jzavala/ETHOS/%s/snapdir_%s/snap_%s.%s.hdf5' % (name,numb,numb,k))
    print file.keys()
    print file['Header'].attrs['Redshift']
    

    print file['Header'].attrs()
    #print file['Header'].attrs['MassTable']
    sys.exit()

    #print list(file['PartType1'].attrs)
    #print list(file['PartType1'].keys())    
    #mass = file['PartType1']['Masses']
    pos = file['PartType1']['Coordinates']

    """file1 = open('/n/dvorkin_lab/anadr/%s_%s_parttype1_%s.txt' % (name,numb,k), 'w')
    
    for i in pos:
        file1.write('%s %s %s\n' % (i[0],i[1],i[2]))
    file1.close()"""




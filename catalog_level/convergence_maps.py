from __future__ import division
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import pylab as py
import sys
import time
from functions import *
from scipy.interpolate import interp1d
import argparse
from tburk_convergence import *
from tnfw_convergence import * 

h = 0.6909

parser = argparse.ArgumentParser()

parser.add_argument('--proj_file',
                  help='path for file with projections, output of make_projections.py',
                  type=str)
parser.add_argument('-o','--outdir',
                    default='./',
                    help='output directory',
                    type=str)
parser.add_argument("-p","--pix_num",
                    default=1011,
                    help="number of pixels in the image",
                    type=int)
parser.add_argument("-s","--side",
                    default=100,
                    help='physical size of the image in kpc',
                    type=int)
parser.add_argument("--name",
                    default='CDM',
                    help='which DM model to use',
                    type=str)
parser.add_argument("-z","--z",
                    default=0.5,
                    help="simulation redshift")
parser.add_argument("--m_high",
                   default=1e-2,
                   help='highest subhalo mass, in units of 10^10 M_sun')
parser.add_argument("--m_low",
                   default=1e-4,
                   help='highest subhalo mass, in units of 10^10 M_sun')
parser.add_argument("-n","--num_proj",
                   default=10,
                   help='total number of projections divided by 3')
parser.add_argument('-S','--sigma_crit',
                    default=2.35e9,
                    help='critical surface mass density for lensing, in units of M_sun/kpc^2',
                    type=float)

args = parser.parse_args()
proj_file = args.proj_file
outdir = args.outdir
pix_num = args.pix_num
rnge = args.side
name = args.name
z = args.z
mhigh = args.m_high
mlow = args.m_low
num_proj = args.num_proj
sigmac = args.sigma_crit

######################################################

if z == 0.5:
    numb = '095'

elif z == 0:
    numb='127'
    z = 0.5 #redshift changed to give physical values  

elif z == 1:
    numb = '078'

else:
    print("wrong redshift specified")
    sys.exit()

######################################################

if mhigh >= 10:
    mlab = 'm11'

elif mhigh == 1:
    mlab = 'm10'

elif mhigh ==1e-1:
    mlab = 'm9'

elif mhigh == 1e-2:    
     mlab = 'm8'

elif mhigh == 1e-3:
    mlab = 'm7'

######################################################

if mlow == 1e-4:
    
    mpart = 1.90382332e-6 / h
    m_min = 50 * mpart
    mlab2 = 'm6'

elif mlow == 1e-3:
    
    mlab = 'm7'
    m_min = mlow

elif mlow == 1e-2:
    
    mlab = 'm8'
    m_min = mlow

elif mlow == 1e-1:
 
    mlab = 'm9'
    m_min = mlow

elif mlow == 1:

    mlab = 'm10'
    m_min = mlow

elif mlow == None:
    mpart = 1.90382332e4 / h
    m_min = 50 * mpart
    mlab2 = 'm6'
    m_min = mlow

if mlab == mlab2 or mhigh <= m_min:
    print 'wrong mass bounds specified!'
    sys.exit()

######################################################

list1 = []
for i in range(16):
    list1.append('/n/hernquistfs3/jzavala/ETHOS/%s/groups_%s/fof_subhalo_tab_%s.%s.hdf5' % (name,numb,numb,i))


posmass = np.load(proj_file)
print 'loaded projections %s' % proj_file

pos_2d = []
masses = []
rmaxs = []
rhs = []
for i in posmass:
    r,m,rh,rmax,_ = zip(*i)
    pos_2d.append(r)
    masses.append(m)
    rmaxs.append(rmax)
    rhs.append(rh)

masses = np.array([np.array(i) for i in masses])
rhs = np.array([np.array(i) for i in rhs])
rmaxs = np.array([np.array(i) for i in rmaxs])

#############################

taus = []
mscale = []

if name == 'CDM':

    tau_vs_xhalf_arr = np.loadtxt('tau_vs_xhalf.txt')
    tau_vs_xhalf = interp1d(tau_vs_xhalf_arr[:,0],tau_vs_xhalf_arr[:,1])

    r_s = rmaxs/2.1626


    for i,j,k in zip(rhs,r_s,masses):
 
        tau = tau_vs_xhalf(i/j)
        m_NFW = k * ((tau**2 + 1)**2/(tau**2*((tau**2 - 1)*np.log(tau) + np.pi*tau - (tau**2 + 1))))

        taus.append(np.array(tau))
        mscale.append(np.array(m_NFW))

elif name == 'ETHOS_4':

    tau_vs_xhalf_arr = np.loadtxt('tau_vs_xhalf_burk.txt')
    tau_vs_xhalf = interp1d(tau_vs_xhalf_arr[:,0],tau_vs_xhalf_arr[:,1])

    r_b = rmaxs/3.3446
    p = 0.66651

    for i,j,k in zip(rhs,r_b,masses):

	select1 = (i/j) < min(tau_vs_xhalf_arr[:,0])
	select2 = (i/j) > max(tau_vs_xhalf_arr[:,0])
	
	tau = np.zeros( (i/j).shape )
 	tau[~(select1 ^ select2)] = tau_vs_xhalf( (i/j)[~(select1 ^ select2)] )
	
	ex = extrap1d(tau_vs_xhalf)
 	tau[(select1 ^ select2)] = ex( (i/j)[(select1 ^ select2)] )

        n = (tau**2 * (np.pi * (p-tau)**2 + 4*tau**2 * np.log(p/tau)))/(4 * (p**4 - tau**4))
        m_burk = k / n
        
	taus.append(np.array(tau))
        mscale.append(np.array(m_burk))

    r_s = r_b / p 

r_t = np.array(taus) * r_s

shiftx = 0
shifty = 0
pix_size = rnge/pix_num
x = np.linspace(-rnge/2+shiftx, rnge/2+shiftx, pix_num)
y = np.linspace(-rnge/2+shifty, rnge/2+shifty, pix_num)

count = 0
conv_list = []

for g,q,u,z in zip(pos_2d,taus,mscale,r_s):

    count += 1
    print count
    
    beginning = time.time()
    
    if name == 'CDM':
        print 'using the truncated NFW profile'
        result = [conv_tnfw(x[:,None]-i[0],y[None,:]-i[1],k,l,m,sigmac).real for i,k,l,m in zip(g,q,u,z)]
    
    elif name == 'ETHOS_4':
	print 'using the truncated Burkert profile'
	result = [conv_tburk(x[:,None]-i[0],y[None,:]-i[1],m,l,p,k,sigmac).real for i,k,l,m in zip(g,q,u,z)]

    res = sum(result)
    conv_list.append(res)
    
    print("%s seconds" % (time.time() - beginning))

    if count % 30 == 0:
        
	conv_list2 = np.asarray(conv_list)
	filename = 'conv%s_%s_%s_%s_%s' % (name,numb,pix_num,count,rnge)
        dir3 = outdir + filename
        np.save(dir3,conv_list2)
        print 'saved convergence maps as %s' % dir3
        
	conv_list = []

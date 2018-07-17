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

#posmass = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/projections_%s_%s_%s.npy' % (name,numb,mlab,name,numb,mlab))
posmass = np.load('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/origproj_%s_%s_%s.npy' % (name,numb,mlab,name,numb,mlab))
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
rmaxs = []
r_3d = []
for i in posmass:
    r,m,rh,rmax,r3d = zip(*i)
    pos_2d.append(r)
    masses.append(m)
    rmaxs.append(rmax)
    r_3d.append(r3d)

def rs(msub):
    return 0.1*(msub/10**6)**(1/3)

def mb(m,tau,p):
    n = np.asarray([(i**2 * (np.pi * (p-i)**2 + 4*i**2 * np.log(p/i)))/(4 * (p**4 - i**4)) for i in tau])
    m_b = m/n
    return m_b

def term1(x,p,tau):
    return (2*p*np.sqrt(1/(tau**2+x**2)))/(p**4 - tau**4)

def term2(x,p,tau):
    down = x < p
    up = np.logical_not(down)
    lowx = np.sqrt(1/(x[up]**2-p**2))/(p*(p**2+tau**2))
    highx = 0
    ans = np.zeros(x.shape)
    ans[up] = lowx
    ans[down] = highx
    return ans

def term3(x,p,tau):
    return np.sqrt(1/(x**2+p**2))/(p**3-p*tau**2)

def term5(x,p,tau):
    return 2*np.arctanh(p/np.sqrt(x**2+p**2)) / ((p**3-p*tau**2)*np.sqrt(x**2+p**2))

def term6(x,p,tau):
    return 4*tau*np.arctanh(tau/np.sqrt(x**2+tau**2)) / ((p**4-tau**4)*np.sqrt(x**2+tau**2))

def term4(x,p,tau):
    down = x < p
    up = np.logical_not(down)
    lowx = 2*np.arctan(p/np.sqrt(x[up]**2-p**2)) / ((p**3+p*tau**2)*np.sqrt(x[up]**2 - p**2))
    highx = np.real(2*np.arctan(p/(1j*np.sqrt(p**2-x[down]**2))) / ((p**3+p*tau**2)*1j*np.sqrt(p**2-x[down]**2)))
    ans = np.zeros(x.shape)
    ans[up] = lowx
    ans[down] = highx
    return ans

def conv_tburk_2d(x,y,rs,mb,pval,tauval):
    p = pval
    t = tauval
    z = np.sqrt(x**2+y**2)/rs
    term_1 = term1(z,p,t)
    term_2 = term2(z,p,t)
    term_3 = term3(z,p,t)
    term_4 = term4(z,p,t)
    term_5 = term5(z,p,t)
    term_6 = term6(z,p,t)
    sig = (mb/(8*np.pi*sigmac*rs**2))*t**2 * (np.pi*(term_1 - term_2 - term_3) + term_4 - term_5 + term_6)
    return sig

masses = np.asarray([np.asarray(i) for i in masses])
r_s = rs(masses)
avg_min_rs = np.mean([min(i) for i in r_s])

def rt_func(m,r3d):
    return (m/10**6)**(1/3) * (r3d/100)**(2/3)

r3d_arr = np.asarray([np.asarray(i) for i in r_3d])
rt = rt_func(masses,r3d_arr)

max_rt = [max(i) for i in rt]
avg_max_rt = np.mean(max_rt)

p = 0.7
t = rt/r_s
m_b = mb(masses,t,p)

rnge = 100
shiftx = 0
shifty = 0
pix_num = 501
pix_size = rnge/pix_num
x = np.linspace(-rnge/2+shiftx, rnge/2+shiftx, pix_num)
y = np.linspace(-rnge/2+shifty, rnge/2+shifty, pix_num)

count = 0
conv_list = []

for g,q,u,z in zip(pos_2d,t,m_b,r_s):
    count += 1
    print count
    beginning = time.time()
    result = [conv_tburk_2d(x[:,None]-i[0],y[None,:]-i[1],m,l,p,k) for i,k,l,m in zip(g,q,u,z)]
    res = sum(result)
    conv_list.append(res)
    print("%s seconds" % (time.time() - beginning))

    if count % 10 == 0:
        conv_list2 = np.asarray(conv_list)
        np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburk%s_%s_%s_%s_%s_%s' % (name,numb,mlab,name,numb,pix_num,count,rnge,mlab),conv_list2)
        conv_list = []

feat = {'pixels':pix_num,'rnge':rnge,'num_proj':len(posmass),'m_high cut':m_high_cut,'m_high':mhigh}

np.save('/n/home04/adiazrivero/conv_powerspec/clean_fof_masscut/%s/%s/%s/convburkfeat%s_%s_%s_%s_%s' % (name,numb,mlab,name,numb,pix_num,rnge,mlab), feat)

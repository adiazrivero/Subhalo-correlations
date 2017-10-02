import numpy as np
import matplotlib as m
import matplotlib.pyplot as plt
import pylab as py
import sys
import time
from coords_50x50_2 import coords_list

start_time = time.time()

"""avg = []
for i in coords_list:
    avg.append(len(i))
print "Average # of subhalos per projection: %s" % np.mean(avg)"""

bin_avg = []
coadd_ps = []
coadd_corr = []

for i in coords_list:
    x,y = zip(*i)
    H,xedges,yedges = np.histogram2d(x,y,bins=20)
    Nbar = H.mean()
    bin_avg.append(Nbar)
    corr = (H - Nbar) / Nbar
    coadd_corr.append(corr)
    ft = np.fft.fft2(corr)
    ft = [i/(2*np.pi) for i in ft] #numpy fft convention: e^(-i2\pixk)
    ps2D = np.abs(ft)**2
    coadd_ps.append(ps2D)
    #kx = np.fft.fftfreq(xedges.size,d=5)
    #ky = np.fft.fftfreq(yedges.size,d=5)

tot_corr = sum(coadd_corr)
tot_corr = [i/len(coadd_corr) for i in tot_corr]

"""py.figure('(n - nbar)/nbar')
plt.imshow(tot_corr,extent=[-50, 50, -50, 50],interpolation='nearest')
plt.colorbar()
py.show()"""

coadd_ps = [i/len(coadd_ps) for i in coadd_ps]
tot_ps = sum(coadd_ps)
tot_ps = np.fft.fftshift(tot_ps)

pix_val = np.arange(-50.1,50.1,5)
pix_val_k = [(2*np.pi)/i for i in pix_val]
pix_2 = [i for i in pix_val_k if i >= 0]

py.figure('2d Power Spectrum')
py.clf()
py.imshow(tot_ps,extent=[pix_val_k[0],pix_val_k[-1], pix_val_k[0],pix_val_k[-1]],interpolation='nearest')
plt.colorbar()
py.show()

"""py.figure('2d Power Spectrum (2)')
py.clf()
py.imshow(tot_ps,extent=[min(kx),max(kx), min(kx),max(kx)],interpolation='nearest')
plt.colorbar()
py.show()"""

def oneD_ps(data,pixel_size=1):
    ps2d = np.array(data)
    pixx, pixy = ps2d.shape
    x1 = np.arange(-pixx/2.,pixx/2.)
    y1 = np.arange(-pixy/2.,pixy/2.)
    x,y = np.meshgrid(y1,x1)

    k = np.sqrt(x**2+y**2)
    kmax = k.max()
    dk = pixel_size
    K = np.arange(kmax/dk)*dk + dk/2.
    nk = len(K)

    ring_mean = []
    for i in range(nk):
        kmin = i*dk
        kmax = kmin + dk
        ind = (k>=kmin) * (k<kmax)
        #print ind*1
        ring_mean.append(data[ind].mean())

    return ring_mean,K

pix_size_k = 2*np.pi*0.2 # 0.2 = 1/(pixel size in r space), which = 5 kpc/h here.
ps1d,K = oneD_ps(tot_ps,pixel_size=pix_size_k)
plt.loglog(K,ps1d)
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

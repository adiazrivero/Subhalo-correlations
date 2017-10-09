import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time

start_time = time.time()

Rnge = 100
num = 1406 # average number of ETHOS subhalos within r < 50 kpc/h after rotating and projecting

coords = []
count = 0
while count < 3000:
    count += 1
    rand = np.random.uniform(-Rnge/2.,Rnge/2.,(num,2))
    coords.append(rand)

rnge = 100
bns = 20
bin_size = rnge/bns
bin_avg = []
coadd_ps = []
coadd_corr = []

for i in coords:
    x,y = zip(*i)
    H,xedges,yedges = np.histogram2d(x,y,bins=bns)
    Nbar = H.mean()
    bin_avg.append(Nbar)
    corr = (H - Nbar) / Nbar
    coadd_corr.append(corr)
    ft = np.fft.fft2(corr)
    ft = [i/(2*np.pi) for i in ft] #numpy fft convention: e^(-i2\pixk)
    ps2D = np.abs(ft)**2
    coadd_ps.append(ps2D)

print "average nbar: %s " % np.mean(bin_avg)

tot_corr = sum(coadd_corr)
tot_corr = [i/len(coadd_corr) for i in tot_corr]

py.figure('(n - nbar)/nbar')
plt.imshow(tot_corr,extent=[-rnge/2, rnge/2, -rnge/2, rnge/2],interpolation='nearest')
plt.colorbar()
py.show()

coadd_ps = [i/len(coadd_ps) for i in coadd_ps]
tot_ps = sum(coadd_ps)
tot_ps = np.fft.fftshift(tot_ps)
kx = np.fft.fftfreq(bns,d=bin_size)

"""py.figure('2d Power Spectrum (2)')
py.imshow(tot_ps,extent=[min(kx),max(kx),min(kx),max(kx)],interpolation='nearest')
plt.colorbar()
py.show()"""

def oneD_ps(data,pixel_size=1):
    ps2d = np.array(data)
    pixx, pixy = ps2d.shape
    x1 = np.arange(-pixx/2.,pixx/2.)
    y1 = np.arange(-pixy/2.,pixy/2.)
    x,y = np.meshgrid(y1,x1)

    conv = (pixx/2.)/kx[1]

    k = np.sqrt(x**2+y**2)
    kmax = k.max()
    dk = pixel_size
    K = (np.arange(kmax/dk)*dk + dk/2.)/conv
    nk = len(K)

    ring_mean = []
    for i in range(nk):
        kmin = i*dk
        kmax = kmin + dk
        ind = (k>=kmin) * (k<kmax)
        #print ind*1
        ring_mean.append(data[ind].mean())

    return ring_mean,K

pix_size_k = 2*np.pi / bin_size
ps1d,K = oneD_ps(tot_ps,pixel_size=pix_size_k)
ps1d_units = [25*i for i in ps1d] # pixel area = 25 (kpc/h)^2
plt.loglog(K,ps1d_units)
plt.xlim(min(K),max(K))
plt.ylim(min(ps1d_units),max(ps1d_units))
plt.show()

print("--- %s seconds ---" % (time.time() - start_time))

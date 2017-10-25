from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import sys
import time

def test(rnge=100,shift=0,subh_num=1,show_nr=False,show_ps=False,num_proj=3000):
    #subh_num = 1406 # average number of ETHOS subhalos within r < 50 kpc/h after rotating and projecting
    #num = 430 # average number of ETHOS subhalos within r < 50 kpc/h after rotating and projecting

    coords = []
    count = 0
    while count < num_proj:
        count += 1
        rand = np.random.uniform(-rnge/2.,rnge/2.,(subh_num,2))
        coords.append(rand)

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
        ft = 1/(2*np.pi)*np.fft.fft2(corr)
        ps2D = np.abs(ft)**2
        coadd_ps.append(ps2D)

    tot_corr = sum(coadd_corr)
    tot_corr = [i/len(coadd_corr) for i in tot_corr]

    if show_nr == True:
        py.figure('(n - nbar)/nbar')
        plt.imshow(tot_corr,extent=[-rnge/2+shift, rnge/2+shift, -rnge/2+shift, rnge/2+shift],interpolation='nearest')
        plt.colorbar()
        py.show()

    coadd_ps = [i/len(coadd_ps) for i in coadd_ps]
    tot_ps = sum(coadd_ps)
    tot_ps = np.fft.fftshift(tot_ps)
    kx = 2*np.pi*np.fft.fftfreq(tot_ps.shape[0],d=bin_size)
    kx = np.fft.fftshift(kx)
    ky = 2*np.pi*np.fft.fftfreq(tot_ps.shape[1],d=bin_size)
    ky = np.fft.fftshift(ky)

    if show_ps == True:
        py.figure('2d Power Spectrum')
        py.imshow(tot_ps,extent=[min(kx),max(kx),min(ky),max(ky)],interpolation='nearest')
        plt.colorbar()
        py.show()

    def oneD_ps(data,pixel_size=1):
        ps2d = np.array(data)
        x,y = np.meshgrid(ky,kx)
        k = np.sqrt(x**2+y**2)
        kmax = k.max()
        dk = pixel_size
        K = np.arange(kmax/dk)*dk

        ring = []
        for i in range(len(K)):
            kmin = i*dk
            kmax = kmin + dk
            index = (k>=kmin) * (k<=kmax)
            ring.append(data[index].mean())

        return ring,K

    pix_size_k = np.abs(kx[0]-kx[1])
    ps1d,K = oneD_ps(tot_ps,pixel_size=pix_size_k)

    A_pix = bin_size**2
    A_box = rnge**2
    norm = A_pix**2/A_box
    ps1d = [norm*i for i in ps1d]

    return ps1d

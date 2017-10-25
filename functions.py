from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import h5py
import sys

def simul_data(Rnge=1):
    list1 = ['fof_subhalo_tab_127.0.hdf5','fof_subhalo_tab_127.1.hdf5','fof_subhalo_tab_127.2.hdf5','fof_subhalo_tab_127.3.hdf5','fof_subhalo_tab_127.4.hdf5','fof_subhalo_tab_127.5.hdf5','fof_subhalo_tab_127.6.hdf5','fof_subhalo_tab_127.7.hdf5','fof_subhalo_tab_127.8.hdf5','fof_subhalo_tab_127.9.hdf5','fof_subhalo_tab_127.10.hdf5','fof_subhalo_tab_127.11.hdf5','fof_subhalo_tab_127.12.hdf5','fof_subhalo_tab_127.13.hdf5','fof_subhalo_tab_127.14.hdf5','fof_subhalo_tab_127.15.hdf5']

    masses = []
    positions1 = []
    vmax1 = []
    for i in list1:
        file = h5py.File(i,'r')
        if len(file['Subhalo'].keys()) != 0:
            masses.append(file['Subhalo']['SubhaloMass'].value) #in 10^10 M_sun/h
            positions1.append(file['Subhalo']['SubhaloPos'].value) # in kpc/h, box coordinates

    masses = [i for j in masses for i in j]
    positions1 = [i for j in positions1 for i in j]
    parent_mass = masses[0]
    parent_pos = positions1[0]
    positions = [i - parent_pos for i in positions1]

    subh = zip(masses[1:],positions[1:])
    subh_cut=[]
    [subh_cut.append(i) for i in subh if i[0] > 1.5e-4 and -Rnge/2. < i[1][0] < Rnge/2. and -Rnge/2. < i[1][1] < Rnge/2. and -Rnge/2. < i[1][2] < Rnge/2.]

    subh_mass,subh_pos = zip(*subh_cut)
    #print "total number of ETHOS subhalos within r < %s kpc/h: %s " % (Rnge/2,len(subh_cut))

    return subh_mass,subh_pos

def rotation(nx,ny,nz,theta):
    R = [[np.cos(theta) + (nx**2)*(1-np.cos(theta)) , nx*ny*(1-np.cos(theta)) - nz*np.sin(theta) , nx*nz*(1-np.cos(theta)) + ny*np.sin(theta)], [nx*ny*(1-np.cos(theta)) + nz*np.sin(theta) , np.cos(theta) + (ny**2)*(1-np.cos(theta)) , ny*nz*(1-np.cos(theta)) - nx*np.sin(theta)], [nz*nx*(1-np.cos(theta)) - ny*np.sin(theta) , nz*ny*(1-np.cos(theta)) + nx*np.sin(theta) , np.cos(theta) + (nz**2)*(1-np.cos(theta))]]

    return R

def projections(positions,rnge=100,shift=0,num_proj=1000):
    coords = []
    count = 0
    while count < num_proj:
        count += 1
        nnx = np.random.uniform(0,10)
        nny = np.random.uniform(0,10)
        nnz = np.random.uniform(0,10)
        theta = np.random.uniform(0,2*np.pi)

        nx = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnx
        ny = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nny
        nz = 1/np.sqrt(nnx**2 + nny**2 + nnz**2) * nnz

        R = rotation(nx,ny,nz,theta)
        rot_pos = [np.dot(R,i) for i in positions]

        proj_xy = [[i[0],i[1]] for i in rot_pos if -rnge/2.+shift < i[0] < rnge/2.+shift and -rnge/2.+shift < i[1] < rnge/2.+shift]
        proj_xz = [[i[0],i[2]] for i in rot_pos if -rnge/2.+shift < i[0] < rnge/2.+shift and -rnge/2.+shift < i[2] < rnge/2.+shift]
        proj_yz = [[i[1],i[2]] for i in rot_pos if -rnge/2.+shift < i[1] < rnge/2.+shift and -rnge/2.+shift < i[2] < rnge/2.+shift]

        coords.append(proj_xy)
        coords.append(proj_xz)
        coords.append(proj_yz)

        tot_num_subh = []
        for i in coords:
            tot_num_subh.append(len(i))
        avg_num_subh = np.mean(tot_num_subh)
        #print "average number of subhalos within r < 50 kpc/h after rotating & projecting: %s" % np.mean(tot_num_subh)

    return coords,avg_num_subh


def twoD_ps(coords,bns=20,rnge=100,shift=0,show_nr=False,show_ps=False):
    """
    takes in projected subhalo positions and returns the 2D power spectrum:
        coords: array of 2d subhalo positions
        bns: number of bins in the horizontal and vertical axes of the 2D power spectrum
        rnge: spatial range of the subhalo positions
        shift: offset from the halo center
        show_nr: whether you want to see xi_ss
        show_ps: whether you want to see the 2d power spectrum
    """
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

    #print "average nbar: %s " % np.mean(bin_avg)

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
    #print kx
    #print [2*np.pi/i for i in kx]

    if show_ps == True:
        py.figure('2d Power Spectrum')
        py.imshow(tot_ps,extent=[min(kx),max(kx),min(ky),max(ky)],interpolation='nearest')
        plt.colorbar()
        py.show()

    return coadd_ps,tot_ps,kx,ky,bin_size

def oneD_ps(data,kx,ky,rnge=100,bin_size=1,pixel_size=1):
    ps2d = np.array(data)
    x,y = np.meshgrid(ky,kx)
    k = np.sqrt(x**2+y**2)
    kmax = k.max()
    dk = pixel_size
    K = np.arange(kmax/dk)*dk
    """print K
    print [2*np.pi/i for i in K]"""

    ps1d = []
    for i in range(len(K)):
        kmin = i*dk
        kmax = kmin + dk
        index = (k>=kmin) * (k<=kmax)
        #print index*1
        ps1d.append(data[index].mean())

    A_pix = bin_size**2
    A_box = rnge**2
    norm = A_pix**2/A_box
    ps1d = [norm*i for i in ps1d]

    return ps1d,K,norm

def error_bars(coadd_ps,ps1d,kx,ky,rnge=100,pix_size=1):
    var = []
    for i in coadd_ps:
        ps,kk,norm = oneD_ps(i,kx,ky,rnge=100,bin_size=1,pixel_size=pix_size)
        ps = [norm*i for i in ps]
        variance = [(i-j)**2 for i,j in zip(ps,ps1d)]
        var.append(variance)

    x = 1/len(var)
    #vari = [x*sum(i) for i in zip(*var)]
    error = x*np.sum(var,0)

    return error

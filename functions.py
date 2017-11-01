from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import h5py
import sys


def simul_data(data,Rnge=1):
    """
    takes in simulation data and returns subhalo masses and positions
    inputs:
        data: simulation data
        Rnge: (2*radial distance) out to which we want to include subhalos
    outputs:
        subh_mass: list of subhalo masses
        subh_pos: list of subhalo positions
    """
    masses = []
    positions1 = []
    vmax1 = []
    for i in data:
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

    return subh_mass,subh_pos

def rotation(nx,ny,nz,theta):
    """
    input:
        unit vector (nx,ny,nz)
        theta: rotation angle
    output:
        3x3 rotation matrix that rotates a 3d vector about a unit vector (nx,ny,nz) by an angle theta
    """
    R = [[np.cos(theta) + (nx**2)*(1-np.cos(theta)) , nx*ny*(1-np.cos(theta)) - nz*np.sin(theta) , nx*nz*(1-np.cos(theta)) + ny*np.sin(theta)], [nx*ny*(1-np.cos(theta)) + nz*np.sin(theta) , np.cos(theta) + (ny**2)*(1-np.cos(theta)) , ny*nz*(1-np.cos(theta)) - nx*np.sin(theta)], [nz*nx*(1-np.cos(theta)) - ny*np.sin(theta) , nz*ny*(1-np.cos(theta)) + nx*np.sin(theta) , np.cos(theta) + (nz**2)*(1-np.cos(theta))]]

    return R

def projections(positions,rnge=100,shift=0,num_proj=1000):
    """
    takes in arrays of 3d positions and returns arrays of 2d positions
    inputs:
        positions: array of 3d positions
        rnge: (2*radial distance) out to which we want to include subhalos (i.e. box size)
        shift: optional shift away from the host center; default is zero
        num_proj: (total number of 2d maps we want)/3 (the division by three is because for each 3d map we project in three different directions, namely xy,xz,yz)
    outputs:
        coords: array of 2d positions
        avg_num_subh: average number of subhalos within r < rnge/2 kpc/h after rotating and projecting
    """
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
    inputs:
        coords: list of lists, where each list is an array of 2d subhalo positions
        bns: number of bins in the horizontal and vertical axes of the 2D power spectrum
        rnge: box size
        shift: offset from the halo center
        show_nr: whether you want to display xi_ss
        show_ps: whether you want to display the 2d power spectrum
    outputs: individual_ps,tot_ps,kx,ky,bin_size
        individual_ps: list of lists, where each list is a 2d power spectrum
        tot_ps: total 2D power spectrum after coadding all 2d power spectra
        ky,ky: fft frequencies
        bin_size: pixel size
    """
    bin_size = rnge/bns
    bin_avg = []
    individual_ps = []
    coadd_corr = []

    for i in coords:
        x,y = zip(*i)
        N,xedges,yedges = np.histogram2d(x,y,bins=bns)
        Nbar = N.mean()
        bin_avg.append(Nbar)
        corr = (N - Nbar) / Nbar
        coadd_corr.append(corr)
        ft = np.fft.fft2(corr)
        ft = 1/(2*np.pi)*np.fft.fft2(corr)
        ps2D = np.abs(ft)**2
        individual_ps.append(ps2D)

    print "average Nbar: %s " % np.mean(bin_avg)
    nbar = Nbar/(bin_size**2)
    print "average nbar: %s " % np.mean(nbar)

    tot_corr = sum(coadd_corr)
    tot_corr = [i/len(coadd_corr) for i in tot_corr]

    if show_nr == True:
        py.figure('(n - nbar)/nbar')
        plt.imshow(tot_corr,extent=[-rnge/2+shift, rnge/2+shift, -rnge/2+shift, rnge/2+shift],interpolation='nearest')
        plt.colorbar()
        py.show()

    individual_ps = [i/len(individual_ps) for i in individual_ps]
    tot_ps = sum(individual_ps)
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

    return individual_ps,tot_ps,kx,ky,bin_size

def angular_average(data,kx,ky,rnge=100,bin_size=1,dk=1):
    """
    takes in a 2d map and returns a 1d, angularly-averaged map
    inputs:
        data: 2d map
        ky,ky: fft frequencies
        rnge: 2d map range
        bin_size:
        pixel_size:

    """
    ps2d = np.array(data)
    x,y = np.meshgrid(ky,kx)
    k = np.sqrt(x**2+y**2)
    kmax = max(kx)

    dk = dk
    K = np.arange(kmax/dk)*dk

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
    ps1d = [norm*i for i in ps1d][1:]
    K = K[1:]
    return ps1d,K,norm

def variance(individual_ps,ps1d,N,kx,ky,rnge=100,pix_size=1):
    """
    takes in a list of lists where each nested list is a single 2d map and returns the variance of each map with respect to the average of all the maps
        individual_ps: list of lists where each nested list is a single 2d map
        ps1d: average of all 2d maps
        ky,ky: fft frequencies
        rnge: spatial extent out to which we want to consider subhalos
        pix_size: pixel size
    """
    var = []
    for i in individual_ps:
        ps,kk,norm = angular_average(i,kx,ky,rnge=100,bin_size=1,dk=pix_size)
        variance = [(i-j)**2 for i,j in zip(ps,ps1d)]
        var.append(variance)
    x = 1/N
    var = x*np.sum(var,0)
    return var

def poisson_realization(rnge=100,subh_num=1,num_proj=3000):
    """
    Creates 2d maps with a given number of randomly distributed subhalos
        rnge: box size
        subh_num: desired number of subhalos in the map
        num_proj: desired number of projections, default is 3000
    """
    coords = []
    count = 0
    while count < num_proj:
        count += 1
        rand = np.random.uniform(-rnge/2.,rnge/2.,(int(subh_num),2))
        coords.append(rand)
    return coords

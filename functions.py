from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pylab as py
import h5py
import sys
import itertools
import scipy
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from fullprint import fullprint

def simul_data(data,Rnge=600,mhigh_cut=False):
    """
    takes in subhalo catalogs and returns subhalo masses, positions, halfmass and vmax radii
    inputs:
        data: halo catalog
        Rnge: size (in kpc/h) of the box/cube out to which we want to include subhalos
    outputs:
        subh_mass: list of subhalo masses
        subh_pos: list of subhalo positions
        rh : halfmass radius
        rvmax: radius at which v_max (max. circular velocity) is achieved
    """
    masses = []
    positions1 = []
    rh = []
    rvmax = []
    for i in data:
        file = h5py.File(i,'r')
        if len(file['Subhalo'].keys()) != 0:
            masses.append(file['Subhalo']['SubhaloMass'].value) #in 10^10 M_sun/h
            positions1.append(file['Subhalo']['SubhaloPos'].value) # in kpc/h, absolute box coordinates
            rh.append(file['Subhalo']['SubhaloHalfmassRad'].value) # in kpc/h
            rvmax.append(file['Subhalo']['SubhaloVmaxRad'].value) # in kpc/h

    masses = [i for j in masses for i in j]
    positions1 = [i for j in positions1 for i in j]
    rh = [i for j in rh for i in j]
    rvmax = [i for j in rvmax for i in j]
    parent_mass = masses[0]
    parent_pos = positions1[0]
    positions = [i - parent_pos for i in positions1] #setting origin at halo center
    subh = zip(masses[1:],positions[1:],rh[1:],rvmax[1:])

    #keeping subhalos above mass resolution limit (and below mhigh if mhigh_cut = True) and out to 300 kpc/h

    if mhigh_cut == True:

        subh_cut = [i for i in subh if 1.5e-4 < i[0] <= 1e-2 and -Rnge/2. < i[1][0] < Rnge/2. and -Rnge/2. < i[1][1] < Rnge/2. and -Rnge/2. < i[1][2] < Rnge/2.] #imposing an mhigh cut
        mass,pos,rh,rvmax = zip(*subh_cut)

        return mass,pos,rh,rvmax

    else:

        subh_cut = [i for i in subh if 1.5e-4 < i[0] and -Rnge/2. < i[1][0] < Rnge/2. and -Rnge/2. < i[1][1] < Rnge/2. and -Rnge/2. < i[1][2] < Rnge/2.]
        mass,pos,rh,rvmax = zip(*subh_cut)

        return mass,pos,rh,rvmax

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

#THIS ONE WORKS WITH COORDSM
def projections(positions,masses,rnge=100,shift=0,num_proj=1000):
    coords = []
    coordsm = []
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

        proj_xy2 = [[[i[0],i[1]],j] for i,j in zip(rot_pos,masses) if -rnge/2.+shift < i[0] < rnge/2.+shift and -rnge/2.+shift < i[1] < rnge/2.+shift]
        proj_xz2 = [[[i[0],i[2]],j] for i,j in zip(rot_pos,masses) if -rnge/2.+shift < i[0] < rnge/2.+shift and -rnge/2.+shift < i[2] < rnge/2.+shift]
        proj_yz2 = [[[i[1],i[2]],j] for i,j in zip(rot_pos,masses) if -rnge/2.+shift < i[1] < rnge/2.+shift and -rnge/2.+shift < i[2] < rnge/2.+shift]

        coordsm.append(proj_xy2)
        coordsm.append(proj_xz2)
        coordsm.append(proj_yz2)

        tot_num_subh = []
        for i in coords:
            tot_num_subh.append(len(i))
        avg_num_subh = np.mean(tot_num_subh)

    return coordsm,avg_num_subh

def twoD_nr(coords,pix_num=20,rnge=100,shift=0,show_nr=False):
    """
    takes in projected subhalo positions and returns the 2D n(r)
    inputs:
        coords: list of lists, where each list is an array of 2d subhalo positions
        pix_num: number of bins in the horizontal and vertical axes of the 2D power spectrum
        rnge: box size
        shift: offset from the halo center
        show_nr: whether you want to display it; default is False
    outputs:
        ind_corr: list of lists, where each list is a 2d n(r) map
        tot_corr: coadded 2d n(r) map (i.e. coadding maps in ind_corr)
        bin_size: pixel size
        np.mean(nbar): average nbar across all maps
    """
    bin_size = rnge/pix_num
    bin_avg = []
    ind_corr = []
    for i in coords:
        x,y = zip(*i)
        n,xedges,yedges = np.histogram2d(x,y,bins=pix_num)
        N = n.T
        Nbar = N.mean()
        bin_avg.append(Nbar)
        corr = (N - Nbar) / Nbar
        ind_corr.append(corr)

    nbar = Nbar/(bin_size**2)

    tot_corr = sum(ind_corr)/len(ind_corr)
    #tot_corr = [i/len(ind_corr) for i in tot_corr]

    if show_nr == True:
        py.figure('(n - nbar)/nbar %s x %s' % (pix_num,pix_num))
        plt.imshow(tot_corr,extent=[-rnge/2+shift, rnge/2+shift, -rnge/2+shift, rnge/2+shift],interpolation='nearest')
        plt.colorbar()
        py.show()

    return ind_corr,tot_corr,np.mean(nbar)

def twoD_ps(data=None,ind_data=None,pix_size=0,rnge=100,shift=0,show_ps=False):
    """
    takes in a 2D array and returns the 2D FFT:
    inputs:
        ind_coords: n 2d arrays, whose average = coadd_coords
        coadd_coords: a singe 2d array
        pix_size: pixel size
        rnge: box size
        shift: offset from the halo center
        show_ps: whether you want to display the 2d power spectrum
    outputs:
        ind_ps_x: list of lists, where each list is a 2d power spectrum
        tot_ps: total 2D power spectrum after coadding all PS in ind_ps_x
        ky,ky: fft frequencies
    """
    if data == None:

        ind_ps = []
        for i in ind_data:
            ft = np.fft.fft2(i)
            ps2D = np.abs(ft)**2
            ind_ps.append(ps2D)

        A_pix = pix_size**2
        A_box = rnge**2
        norm = A_pix**2/A_box

        ind_ps_x = [norm*np.fft.fftshift(i) for i in ind_ps]
        tot_ps = sum(ind_ps_x)/len(ind_ps_x)

        kx = 2*np.pi*np.fft.fftfreq(tot_ps.shape[0],d=pix_size)
        kx = np.fft.fftshift(kx)
        ky = 2*np.pi*np.fft.fftfreq(tot_ps.shape[1],d=pix_size)
        ky = np.fft.fftshift(ky)

        if show_ps == True:
            py.figure('2d Power Spectrum')
            py.imshow(tot_ps,extent=[min(kx),max(kx),min(ky),max(ky)],interpolation='nearest')
            plt.colorbar()
            py.show()

        return ind_ps_x,tot_ps,kx,ky

    elif ind_data == None:

        ft = np.fft.fft2(data)
        ps2D = np.abs(ft)**2
        tot_ps = np.fft.fftshift(ps2D)

        A_pix = pix_size**2
        A_box = rnge**2
        norm = A_pix**2/A_box

        tot_ps = np.asarray([norm*i for i in tot_ps])

        kx = 2*np.pi*np.fft.fftfreq(tot_ps.shape[0],d=pix_size)
        kx = np.fft.fftshift(kx)
        ky = 2*np.pi*np.fft.fftfreq(tot_ps.shape[1],d=pix_size)
        ky = np.fft.fftshift(ky)

        if show_ps == True:
            py.figure('2d Power Spectrum')
            py.imshow(tot_ps,extent=[min(kx),max(kx),min(ky),max(ky)],interpolation='nearest')
            plt.colorbar()
            py.show()

        return tot_ps,kx,ky

def angular_average(data,x,y,rnge=100,pix_num=21,dr=1,remove_first=False):
    """
    takes in a 2d map and returns a 1d, angular-averaged map
    inputs:
        data: 2d map
        ky,ky: fft frequencies
        rnge: 2d map range in kpc/h
        dk: pixel size
    """
    X,Y = np.meshgrid(x,y)
    r = np.sqrt(X**2+Y**2)

    if max(x) <= max(y):
        rmax = max(x)
    else:
        rmax = max(y)
    R = np.arange(rmax/dr)*dr

    #loading a mask that allows us to do the angular averaging
    orig_mask = np.load('mask.npy')
    mask = []
    #keep only the right number of rings (depends on pixel number)
    for i in orig_mask[:len(R)]:
        mask2 = []
        beg = int((len(i[0])-pix_num)/2)
        end = int(len(i[0])-beg)
        for j in i[beg:end]:
            mask2.append(j[beg:end])
        mask.append(mask2)

    data = np.asarray(data)
    ps1d = []
    for i,j in zip(range(len(R)),mask):
        ring = data*j
        ring = ring[ring != 0]
        ps1d.append(np.asarray(ring).mean())

    if remove_first == True:
        return ps1d[1:],R[1:]
    else:
        return ps1d,R

def variance(individual_ps,ps1d,N,kx,ky,rnge=100,pix_num=21,pix_size=1):
    """
    takes in a list of lists where each nested list is a single 2d map and returns the variance of each map with respect to the average of all the maps
        individual_ps: list of lists where each nested list is a single 2d map
        ps1d: average of all 2d maps
        kx,ky: fft frequencies
        rnge: spatial extent out to which we want to consider subhalos
        pix_num: number of pixels
        pix_size: pixel size
    """
    var = []
    for i in individual_ps:
        ps,kk = angular_average(i,kx,ky,rnge=rnge,pix_num=pix_num,dr=pix_size,remove_first=True)
        diff = [(1/N)*(i-j)**2 for i,j in zip(ps,ps1d)]
        var.append(diff)

    var = np.sum(var,axis=0)
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

def interpolation(data,x,y,num=100j,display=False):
    """
    Takes in discrete 2d array and returns an interpolated array
    """
    X, Y = np.meshgrid(x,y)
    l1 = [i for j in X for i in j]
    l2 = [i for j in Y for i in j]
    points = [list(i) for i in zip(l1,l2)]
    values = np.array([i for j in data for i in j])
    #grid_x,grid_y = np.mgrid[-0.5:0.5:num, -0.5:0.5:num]
    grid_x,grid_y = np.mgrid[x[0]:x[-1]:num, y[0]:y[-1]:num]
    grid = griddata(points,values,(grid_x,grid_y),method='cubic')

    if display == True:
        py.figure('Interpolated map')
        plt.imshow(grid.T)
        plt.colorbar()
        plt.show()

    return grid

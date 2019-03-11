Most of the heavy duty functions are in functions.py, tburk_convergence.py and tnfw_convergence.py. Before running anything else, run make_mask.py: will return a npy file called mask_x.npy where x = pix_num, and it is a mask needed to do the azimuthal averaging. The default pix_num is 1011 - to specify a different number of pixels per side run make_mask.py pix_num=XXX, where XXX is the desired number of pixels (int or float).

Then,

1. make_projections.py: returns a npy file called projections_%s_%s_%s_%s.npy % (name,numb,mlab2,rnge), which contains all the different line of sight projections. Arguments that must be passed when compiling (python make_projections.py name=XXX numb=YYY etc.)

    name # if 0 --> CDM simulation // if 2 --> ETHOS_4 simulation
    numb # if 78 --> redshift=1 // if 95 --> redshift=0.5 // if 127 --> redshift=0
    rnge # distance (in kpc) out to which you want to keep subhalos
    mhigh # highest subhalo mass, in units of 1e10 M_sun
    mlow # lowest subhalo mass, in units of 1e10 M_sun --> mlab2 in the file name corresponds to the lowest subhalo mass, e.g. if the lowest subhalo mass is 10^6 M_sun, mlab2=m6. 
    num_proj #number of projections to be done; actual number is num_proj*3 (after each rotation the simulation is projected along the x, y and z axes)

2. convergence_maps.py: makes a tnfw/tburk convergence map for each projection in projections.npy. Arguments that must be passed when compiling (python make_projections.py name=XXX numb=YYY etc.)

    name # if 0 --> CDM simulation // if 2 --> ETHOS_4 simulation
    numb # if 78 --> redshift=1 // if 95 --> redshift=0.5 // if 127 --> redshift=0
    rnge # distance (in kpc) out to which you want to keep subhalos
    mhigh # highest subhalo mass, in units of 1e10 M_sun
    mlow # lowest subhalo mass
    rein # the Einstein radius of the lens
    sigmac # the critical surface mass density for lensing
    zs # redshift of the source
    zl # redshift of the lens
    pix_num # number of pixels
   
Returns them as conv%s_%s_%s_%s_%s_%s.npy %s (name,numb,pix_num,count,rnge,mlab2) count=the count in the loop (they are stored in sets of ten, so for the first 10 maps count, for the next ten count=2, etc.). 

3. power_spectrum.py: it takes in conv_XXX.npy and returns the 1d power spectrum for each individual convergence maps in a file called ind_curves_XXX.npy (where the XXX is as above).

    name # if 0 --> CDM simulation // if 2 --> ETHOS_4 simulation
    numb # if 78 --> redshift=1 // if 95 --> redshift=0.5 // if 127 --> redshift=0
    rnge # distance (in kpc) out to which you want to keep subhalos
    mhigh # highest subhalo mass, in units of 1e10 M_sun
    mlow # lowest subhalo mass
    pix_num # number of pixels per side in the convergence maps
 
4. Plotting functions: code used to make the figures in the paper. Takes in ind_curves_XXX.npy. 

    figure_2.py 
    mass_dependance.py --> Figure 3 bottom panels
    psub_times_sigmac.py // redshift_dependance.py --> Figure 2 top panels (the former for the PRD version and the latter for the arxiv version)
 
5. The directory "host_features" contains text files with features for the host for the CDM and ETHOS_4 simulations at z=0 and z=0.5. The file host_posmass.py can make these files for any desired simulation/redshift. See below for what variables have to be specified when compiling. The directory "rein_sigmac_all" is the same but for the einstein radius and sigma_crit of a given simulation. The file R_ein.py computes them.

---------

An example of a sequence from beginning to end to obtain the 1d power spectrum, for the CDM simulation at z = 0.5 and with a high mass cut of 10^8 M_sun:

python make_projections.py name=0 numb=95 mhigh=10 mlow=1e-4 rnge=100
python convergence_maps.py name=0 numb=95 mhigh=1e-2 mlow=1e-4 rnge=100
python power_spectrum.py name=0 numb=95 mhigh=1e-2 mlow=1e-4 rnge=200

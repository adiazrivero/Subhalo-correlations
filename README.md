1. Run make_mask.py : will return a npy file called mask_x.npy where x = pix_num (default 1011), and it is a mask needed to do the azimuthal averaging

2. Run make_projections.py: returns a npy file called projections.npy which contains all the different line of sight projections. 

3. Run make_convergence_maps.py or make_convergence_maps_tburk.py: makes a tnfw/tburk convergence map for each projection in projections.npy. Returns them in npy files labelled as conv_x_y_z.npy or convburk_x_y_z.npy where x = pix num, z = range, and y = the count in the loop (they are stored in sets of ten, so the first 10 maps are in conv_x_10_z.npy, the next ten in conv_x_20_z.npy, etc.). This also returns convfeat_x_z.npy or convburkfeat_x_z.npy, which has important features of the convergence map.

4. Run projections_features.py: this will print out many features of the projections, basically those corresponding to the tables in the PDF document and will store these in a npy file called projfeat_x_y.npy, where x = pix num and y = image range.

5. Run plot_power_spectrum.py: it takes in convfeat_x_z.npy, projfeat_x_z.npy  and conv_x_y_z.npy and returns the 1d power spectrum. Most of the heavy-duty functions called are in functions.py. If desired the power spectrum and its error bars can be exported as text files. 

6. Extras: 
    1. dNdm.py: plots average subhalo MF + best fit line
    2. n_r.py : plots N(r_3d), N(r_2d) and n(r_2d) *NOTE: THIS CODE IS A MESS FOR NOW*
    
7. The directory "host_features" contains text files with features for the host for the CDM and ETHOS_4 simulations at z=0 and z=0.5. The file host_posmass.py can make these files for any desired simulation/redshift. See below for what variables have to be specified when compiling. The directory "rein_sigmac_all" is the same but for the einstein radius and sigma_crit of a given simulation. The file R_ein.py computes them.


---------

Note that when you run each file you have to specify certain variables:

    1. "name": 
        for CDM name = 0
        for ETHOS_1 name = 1
        for ETHOS_4 name = 2
        
    2. "numb":
        for z = 0 numb = 127
        for z = 0.5 numb = 95
        for z = 1 numb = 78
        
    3. m_high_cut: (bool)
        if True: some high mass threshold is imposed
        if False: no mass cut
        
    4. mhigh: in units of 10^10 M_sun

An example of a sequence from beginning to end, for the CDM simulation at z = 0.5 and with a high mass cut of 10^8 M_sun:

python make_projections.py numb=95 name=0 m_high_cut=True mhigh=1e-2

python make_convergence_maps.py numb=95 name=0 m_high_cut=True mhigh=1e-2

python projections_features.py numb=95 name=0 m_high_cut=True mhigh=1e-2

python plot_power_spectrum.py numb=95 name=0 plot_errors=True save_files=True mhigh=1e-2

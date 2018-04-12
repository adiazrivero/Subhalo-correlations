1. Run make_mask.py : will return a npy file called mask_x.npy where x = pix_num (default 1011), and it is a mask needed to do the azimuthal averaging
2. Run make_projections.py: returns a npy file called projections.npy which contains all the different line of sight projections.  
3. Run make_convergence_maps.py: makes a tnfw convergence map for each projection in projections.npy. Returns them in npy files labelled as conv_x_y_z.npy where x = pix num, z = range, and y = the count in the loop (they are stored in sets of ten, so the first 10 maps are in conv_x_10_z.npy, the next ten in conv_x_20_z.npy, etc.). This also returns convfeat_x_z.npy, which has important features of the convergence map.
4. Run projections_features.py: this will print out many features of the projections, basically those corresponding to the tables in the PDF document and will store these in a npy file called projfeat_x_y.npy, where x = pix num and y = image range.
5. Run plot_power_spectrum.py: it takes in convfeat_x_z.npy, projfeat_x_z.npy  and conv_x_y_z.npy and returns the 1d power spectrum. Most of the heavy-duty functions called are in functions.py. If desired the power spectrum and its error bars can be exported as text files. 
6. Extras: 
    1. dNdm.py: plots average subhalo MF + best fit line
    2. n_r.py : plots N(r_3d), N(r_2d) and n(r_2d) *NOTE: THIS CODE IS A MESS FOR NOW*

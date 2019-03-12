For now I have included the code to reproduce the results in https://arxiv.org/pdf/1809.00004.pdf using the subhalo catalogs (i.e. not at the particle level). 

Most of the heavy duty functions are in functions.py, tburk_convergence.py and tnfw_convergence.py. Before running anything else, run make_mask.py: will return a npy file called mask_x.npy where x = pix_num, and it is a mask needed to do the azimuthal averaging. The default pix_num is 1011 and default size of the image (in kpc) is 100 - to specify different values run it as 

    make_mask.py --pix_num=xxx --side=yyy

Then,

1. make_projections.py: returns a npy file called projections_%s_%s_%s.npy % (name,numb,rnge), which contains all the different line of sight projections. Optional arguments that can be passed when compiling:

        parser.add_argument('-o','--outdir',
                            default='./',
                            help='output directory',
                            type=str)
        parser.add_argument("-p","--pix_num",
                            default=1011,
                            help="number of pixels in the image",
                            type=int)
        parser.add_argument("-s","--side",
                            default=100,
                            help='physical size of the image in kpc',
                            type=int)
        parser.add_argument("--name",
                            default='CDM',
                            help='which DM model to use',
                            type=str)
        parser.add_argument("-z","--z",
                            default=0.5,
                            help="simulation redshift")
        parser.add_argument("--m_high",
                           default=1e-2,
                           help='highest subhalo mass')
        parser.add_argument("--m_low",
                           default=1e-4,
                           help='highest subhalo mass')
        parser.add_argument("-n","--num_proj",
                           default=10,
                           help='total number of projections divided by 3')

2. convergence_maps.py: makes a tnfw/tburk convergence map for each projection in projections.npy. Arguments that must be passed when compiling (python make_projections.py name=XXX numb=YYY etc.)

        parser.add_argument('--proj_file',
                          help='path for file with projections, output of make_projections.py',
                          type=str)
        parser.add_argument('-o','--outdir',
                            default='./',
                            help='output directory',
                            type=str)
        parser.add_argument("-p","--pix_num",
                            default=1011,
                            help="number of pixels in the image",
                            type=int)
        parser.add_argument("-s","--side",
                            default=100,
                            help='physical size of the image in kpc',
                            type=int)
        parser.add_argument("--name",
                            default='CDM',
                            help='which DM model to use',
                            type=str)
        parser.add_argument("-z","--z",
                            default=0.5,
                            help="simulation redshift")
        parser.add_argument("--m_high",
                           default=1e-2,
                           help='highest subhalo mass, in units of 10^10 M_sun')
        parser.add_argument("--m_low",
                           default=1e-4,
                           help='highest subhalo mass, in units of 10^10 M_sun')
        parser.add_argument("-n","--num_proj",
                           default=10,
                           help='total number of projections divided by 3')
        parser.add_argument('-S','--sigma_crit',
                            default=2.35e9,
                            help='critical surface mass density for lensing, in units of M_sun/kpc^2',
                            type=float)
   
Returns them as conv%s_%s_%s_%s_%s.npy %s (name,numb,pix_num,count,rnge) count=the count in the loop (they are stored in sets of ten, so for the first 10 maps count, for the next ten count=2, etc.) in the directory specified in --outdir.  

3. power_spectrum.py: it takes in the convergence maps and returns the 1d power spectrum for each individual convergence map:

        parser.add_argument('--conv_file1',
                          help='path for file with convergence maps, output of convergence_maps.py',
                          type=str)
        parser.add_argument('--conv_file2',
                          help='path for file with convergence maps, output of convergence_maps.py',
                          type=str)
        parser.add_argument('--conv_file3',
                          help='path for file with convergence maps, output of convergence_maps.py',
                          type=str)
        parser.add_argument('--kdir',
                            default='./',
                            help='output directory for wavenumbers',
                            type=str)
        parser.add_argument('--psdir',
                            default='./',
                            help='output directory for the power spectra',
                            type=str)
        parser.add_argument("-p","--pix_num",
                            default=1011,
                            help="number of pixels in the image",
                            type=int)
        parser.add_argument("-s","--side",
                            default=100,
                            help='physical size of the image in kpc',
                            type=int)
        parser.add_argument("--name",
                            default='CDM',
                            help='which DM model to use',
                            type=str)
        parser.add_argument("-z","--z",
                            default=0.5,
                            help="simulation redshift")
        parser.add_argument("--m_high",
                           default=1e-2,
                           help='highest subhalo mass, in units of 10^10 M_sun')
        parser.add_argument("--m_low",
                           default=1e-4,
                           help='highest subhalo mass, in units of 10^10 M_sun')
        parser.add_argument("-n","--num_proj",
                           default=10,
                           help='total number of projections divided by 3')
        parser.add_argument('-S','--sigma_crit',
                            default=2.35e9,
                            help='critical surface mass density for lensing, in units of M_sun/kpc^2',
                            type=float)
 
returns the power spectra for the individual convergence maps as ind_curves_%s_%s_%s_%s_%s.npy % (name,numb,pix_num,rnge) in psdir together with the corresponding wavenumbers k%s_%s_%s_%s_%s.txt' % (name,numb,pix_num,rnge) in kdir. 
 
 
4. Plotting functions: code used to make the figures in the paper.

    figure_2.py 
    
    mass_dependance.py --> Figure 3 bottom panels
    
    psub_times_sigmac.py // redshift_dependance.py --> Figure 2 top panels (the former for the PRD version and the latter for the arxiv version)
 
5. The directory "host_features" contains text files with features for the host for the CDM and ETHOS_4 simulations at z=0 and z=0.5. The file host_posmass.py can make these files for any desired simulation/redshift. See below for what variables have to be specified when compiling. The directory "rein_sigmac_all" is the same but for the einstein radius and sigma_crit of a given simulation. The file R_ein.py computes them.

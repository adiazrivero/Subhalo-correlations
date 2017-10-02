1. run ethos.py: will return four different files called coords_#x#.py where # is a number that gives the size of the image in kpc/h and is either 10,20,30,50 at the moment.

2. fft.py: imports the file output by ethos.py and outputs the subhalo-subhalo power spectrum - right now the x axis in the power spectrum does not represent the correct k values, instead it represents the distance in pixels from the center. I'll change that soon. 

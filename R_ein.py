from matplotlib import pyplot as plt
import sys, platform, os
import camb
import numpy as np
import pylab as pl
import scipy.stats

pars = camb.CAMBparams()
pars.set_cosmology(H0=70, ombh2=0.02222, omch2=0.1197, mnu=0.06, omk=0, tau=0.078)
pars.InitPower.set_params(ns=0.9655, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0);
results = camb.get_results(pars)

G = 4.43*10**(-15) # pc^3 Msun^-1 yr^-2
c = 0.307 # pc yr^-1

def rein_sigmac(M,zs,zl):
    Dos = results.angular_diameter_distance2(0,zs)
    Dol = results.angular_diameter_distance2(0,zl)
    Dls = results.angular_diameter_distance2(zl,zs)

    DS = Dos * 10**6
    DL = Dol * 10**6
    DLS = Dls * 10**6

    sigmac = (c**2 * DS)/(4*np.pi*G*DL*DLS)
    sigmac_kpc = sigmac*10**6

    rein = np.sqrt((4*G*M*DLS)/(c**2 * DS * DL)) #in rad
    rein_as = rein / (4.848137e-6) #in as
    rein_kpc = rein*DL*10**(-3) # in kpc

    return rein_kpc,sigmac_kpc

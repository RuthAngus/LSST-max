from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from teff_bv import teff2bv_orig
import astropy.constants as const

# Global parameters
Gmin, Gmax = 5200, 6000
Gteff = .5 * (Gmin + Gmax)
Kmin, Kmax = 3700, 5200
Kteff = .5 * (Kmin + Kmax)
Mmin, Mmax = 2400, 3700
# Mmin, Mmax = 2400, 3600
Mteff = .5 * (Mmin + Mmax)
Glogg, Klogg, Mlogg = 4.4, 4.4, 4.4

def M2Teff(m):
    """
    Use a simple scaling law to convert mass to teff.
    """
    return (m/const.Msun)**(2.5/4.) * const.Tsun

def Teff2M(t):
    """
    Use a simple scaling law to convert teff to mass.
    """
    return (t/const.Tsun)**(4./2.5) * const.Msun

def Kroupa_IMF(m, dm, N):
    """
    IMF from Kroupa (2001), arXiv: 0009005
    m: the mass (solar masses).
    dm: the mass interval.
    N: a reference number.
    Returns the number of single stars in the mass interval m to m + dm.
    """
    if .01 < m and m < .08:
        a = .3
    elif .08 < m and m < .5:
        a = 1.3
    elif .5 < m and m < 1.:
        a = 2.3
    elif 1. < m:
        a = 2.3
    return N *

def IMF(t, dt, N):
    """
    IMF from Kroupa (2001), arXiv: 0009005, converted to Teff
    t: the Teff (K)
    dt: the temp interval.
    Returns the number of single stars in the temp interval t to t + dt.
    """

def V2r(V, BV):
    return V - .42*BV + .11

def format(fname, age, feh, nG, nK, nM, mvG, mvK, mvM, e_bv):

    ntot = nG + nK + nM

    # create arrays of temperature and logg, G : K : M
    tf = lambda n, Xmin, Xmax: np.random.uniform(Xmin, Xmax, n)
    teffs = np.concatenate((tf(nG, Gmin, Gmax), tf(nK, Kmin, Kmax),
                            tf(nM, Mmin, Mmax)))
    lg = lambda n, logg: np.ones(n) * logg
    loggs = np.concatenate((lg(nG, Glogg), lg(nK, Klogg), lg(nM, Mlogg)))

    # create arrays of feh, b-v, log(age) and log(teff)
    fehs = np.ones(ntot) * feh
    bvs = teff2bv_orig(teffs, fehs, loggs) - e_bv
    logAges = np.ones(ntot) * np.log10(age)
    logTeff = np.log10(teffs)

    rmag = np.concatenate((V2r(mvG, bvs[:nG]), V2r(mvK, bvs[nG:nG+nK]),
                           V2r(mvM, bvs[nG+nK:])))

    data = np.vstack((logAges, bvs, logTeff, rmag))
    np.savetxt(fname, data.T)
    return logAges, bvs, logTeff, rmag

if __name__ == "__main__":

    # IC 4651
#     format("IC4651", 1.778, -.128, 75, 215, 1500, 14.87, 16.01, 18.91)
    format("IC4651", 1.778, -.128, 5000, 5000, 5000, 14.87, 16.01, 18.91)

    # NGC 5316
#     format("NGC5316", .170, .045, 90, 240, 1650, 16.11, 17.25, 20.15)
    format("NGC5316", .170, .045, 5000, 5000, 5000, 16.11, 17.25, 20.15)

    # NGC 2477
#     format("NGC2477", .822, -.192, 90, 240, 1650, 16.44, 17.58, 20.48)
    format("NGC2477", .822, -.192, 5000, 5000, 5000, 16.44, 17.58, 20.48)

# NGC 5316: 90 G dwarfs, 240 K dwarfs, 1650 M dwarfs
# Age = .170 Gyrs, distance = 1208 pc, dm = 10.410
# NGC 2477: 90 G dwarfs, 240 K dwarfs, 1650 M dwarfs
# Age = .822 Gyrs, distance = 1450, dm = 10.807
# IC 4651: 75 G dwarfs, 215 K dwarfs, 1500 M dwarfs.
# Age = 1.778 Gyrs, distance = 888, dm = 9.742

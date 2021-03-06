# coding: utf-8
# # Recovering rotation periods in simulated LSST data

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
from toy_simulator import simulate_LSST
import simple_gyro as sg
import pandas as pd
import sys


def find_nearest(array, value):
    """
    Match a period to a bin.
    array: array of bin heights.
    value: the period of the star.
    Returns the value and index of the bin.
    """
    m = np.abs(array-value) == np.abs(array-value).min()
    return array[m], m


def assign_amps(ps, log10P, log10R, stdR):
    """
    Take periods and bin values and return an array of amplitudes.
    """
    npi = np.array([find_nearest(10**log10P, p) for p in ps])
    nearest_ps, inds = npi[:, 0], npi[:, 1]
    log_ranges = np.array([log10R[i] for i in inds])[:, 0]
    std_ranges = np.array([stdR[i] for i in inds])[:, 0]
    return np.random.randn(len(ps))*std_ranges + log_ranges


def make_arrays(data, temp_bin, ps, teff, rmag):
    """
    Amplitude arrays for each temperature bin
    """
    P, R, std = np.array(data["log10P"]), np.array(data["log10R"]), \
        np.array(data["stdR"])
    if temp_bin == 3500:
        m = teff < 3750
    elif temp_bin == 6000:
        m = teff > 6000
    else:
        m = (temp_bin - 250 < teff) * (teff < temp_bin + 250)
    periods, teffs, rmags = ps[m], teff[m], rmag[m]
    amplitudes = assign_amps(periods, P, R, std)
    return periods, amplitudes, teffs, rmags


def LSST_sig(m):
    """
    Approximate the noise in figure 2 of arxiv:1603.06638 from the apparent
    r-mag.
    Returns the noise in magnitudes and ppm.
    """
    if m < 19:
        return .005
    mags = np.array([19, 20, 21, 22, 23, 24, 25])
    sigs = np.array([.005, .007, .01, .02, .03, .1, .2])
    return sigs[np.abs(mags - m) == np.abs(mags-m).min()][0]

def pgram(N, years, fname):
    ps = np.linspace(2, 100, 1000)  # the period array (in days)

    print("Computing periodograms")
    # Now compute LS pgrams for a set of LSST light curves & save highest peak
    ids = np.arange(N)
    periods = np.zeros_like(ids)
    for i, id in enumerate(ids):
        sid = str(int(id)).zfill(4)
        x, y, yerr = np.genfromtxt("simulations/{0}/{1}.txt".format(fname,
                                   sid)).T
        m = x < years * 365.25
        xt, yt, yerrt = x[m], y[m], yerr[m][m]
        model = LombScargle().fit(xt, yt, yerrt)  # compute pgram
        pgram = model.periodogram(ps)

        # find peaks
        peaks = np.array([j for  j in range(1, len(ps)-1) if pgram[j-1]
                          < pgram[j] and pgram[j+1] < pgram[j]])
        if len(peaks):
            period = ps[pgram == max(pgram[peaks])][0]
        else:
            period = 0

        periods[i] = period
        np.savetxt("results/{0}/{1}_{2}yr_result.txt".format(fname, sid,
                   years), [period])
    np.savetxt("{0}_{1}yr_results.txt".format(fname, years), periods.T)
    return periods

def inject(fname):
    """
    Simulate rotation periods for LSST targets and attempt to recover those
    rotation periods.
    Saves an array of injected periods (days), recovered periods (days), Teff,
    rmag, injected amplitudes (ppm) and noise (ppm).
    'true_ps, periods, logamps, teffs, rmags, true_as, noises_ppm'
    """

    print("Loading Cluster file...")
    # Randomly select targets from a TRILEGAL output.
    logAges, bvs, logTeff, rmag = np.genfromtxt("{0}.dat".format(fname)).T
    teff = 10**logTeff

    # Calculate periods from ages and colours for cool stars
    m = bvs > .4  # select only cool stars
    cool_ages = 10**logAges[m] * 1e-9
    cool_ps = sg.period(cool_ages, bvs[m])
    cool_teffs = teff[m]
    cool_rmags = rmag[m]

    # Draw from a sum of two Gaussians (modelled in another notebook) that
    # describes the period distribution for hot stars. Approximations:
    # I have lumped all stars with colour < 0.4 in together AND I actually
    # used teff = 6250, not B-V = 0.4 in the other notebook.
    hot_ages = 10**logAges[~m] * 1e-9  # select hot stars
    hot_teffs = teff[~m]
    hot_rmags = rmag[~m]

    # copy parameters for two Gaussians from hot_stars ipython notebook
    A1, A2, mu1, mu2, sig1, sig2 = 254.11651209, 49.8149765, 3.00751724, \
        3.73399554, 2.26525979, 8.31739725

    hot_ps = np.zeros_like(hot_ages)
    hot_ps1 = np.random.randn(int(len(hot_ages)*(1 - A2/A1)))*sig1 + mu1
    hot_ps2 = np.random.randn(int(len(hot_ages)*(A2/A1)))*sig2 + mu2
    hot_ps[:len(hot_ps1)] = hot_ps1
    hot_ps[len(hot_ps1):len(hot_ps1) + len(hot_ps2)] = hot_ps2
    tot = len(hot_ps1) + len(hot_ps2)
    hot_ps[tot:] = np.random.randn(len(hot_ps)-tot)*sig2 + mu2

    # combine the modes
    age = np.concatenate((cool_ages, hot_ages))
    ps = np.concatenate((cool_ps, hot_ps))
    teff = np.concatenate((cool_teffs, hot_teffs))
    rmag = np.concatenate((cool_rmags, hot_rmags))

    print("Calculating amplitudes...")
    # Use Derek's results to calculate amplitudes
    # Column headings: log10P, log10R, stdR, Nbin
    d35 = pd.read_csv("data/rot_v_act3500.txt")
    d40 = pd.read_csv("data/rot_v_act4000.txt")
    d45 = pd.read_csv("data/rot_v_act4500.txt")
    d50 = pd.read_csv("data/rot_v_act5000.txt")
    d55 = pd.read_csv("data/rot_v_act5500.txt")
    d60 = pd.read_csv("data/rot_v_act6000.txt")

    # Assign amplitudes
    pers, logamps, teffs, rmags = \
        np.concatenate((make_arrays(d35, 3500, ps, teff, rmag),
                        make_arrays(d40, 4000, ps, teff, rmag),
                        make_arrays(d45, 4500, ps, teff, rmag),
                        make_arrays(d50, 5000, ps, teff, rmag),
                        make_arrays(d55, 5500, ps, teff, rmag)),
                       axis=1)
#                         make_arrays(d60, 6000, ps, teff, rmag)),
    amps = 10**logamps  # parts per million
    noises_mag = np.array([LSST_sig(mag) for mag in rmags])
    noises_ppm = (1 - 10**(-noises_mag/2.5)) * 1e6

    # Simulate light curves
    print("Simulating light curves...")
    path = "simulations/{0}".format(fname)  # where to save the lcs
    [simulate_LSST(i, pers[i], amps[i], path, noises_ppm[i]) for i in
     range(len(pers))]

    # save the true values
    ids = np.arange(len(pers))
    data = np.vstack((ids, pers, amps))
    np.savetxt("{0}/truth.txt".format(path), data.T)

    print("Saving results")
    data = np.vstack((pers, amps, teffs, rmags, noises_ppm))
    np.savetxt("parameters_{0}.txt".format(fname), data.T)
    return pers, amps, teffs, rmags, noises_ppm

if __name__ == "__main__":
    fname = "{0}".format(sys.argv[1])

    # Run simlations
    pers, amps, teffs, rmags, noises_ppm = inject("{0}".format(fname))

    # recover periods
    pers, amps, teffs, rmags, noises_ppm = \
            np.genfromtxt("parameters_{0}.txt".format(fname)).T
    N = len(pers)
    years = [1, 5, 10]
    for year in years:
        periods = pgram(N, year, fname)
        data = np.vstack((pers, periods, np.log(amps), teffs, rmags, amps,
                          noises_ppm))
        np.savetxt("{0}yr_results{1}.txt".format(year, fname), data.T)

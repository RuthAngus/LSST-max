# coding: utf-8
# # Recovering rotation periods in simulated LSST data

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
from toy_simulator import simulate_LSST
import sys
from LSST_inject_and_recover import LSST_sig

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

def simulate_sin(id, b, p, a, noise, path):
    """
    id: id of the star.
    p: rotation period in seconds
    a: amplitude in ppm
    path: path to save files
    """
    id = str(int(id)).zfill(4)

    # The time array
    x, depth = np.genfromtxt("b{0}_cadence.txt".format(b)).T
    noise_free_y = a * np.sin(2*np.pi*(1./p))
    y = noise_free_y + noise*1e-6 * np.random.randn(len(x))
    yerr = np.ones_like(y) * noise * 1e-6

    data = np.vstack((x, y, yerr))
    np.savetxt("{0}/{1}/{2}.txt".format(path, b, id), data.T)
    truths = np.array([p, a])
    np.savetxt("{0}/{1}/{2}_truth.txt".format(path, b, id), truths)

    plt.clf()
    plt.errorbar(x/365.25, y, yerr=yerr, fmt="k.", capsize=0, ecolor=".5")
    plt.xlabel("$\mathrm{Time~(years)}$")
    plt.ylabel("$\mathrm{Normalised~Flux}$")
    plt.xlim(min(x/365.25), max(x/365.25))
    plt.subplots_adjust(left=.2, bottom=.12)
    plt.savefig("test.pdf".format(id))

def inject(fname, N, b):
    """
    Simulate sinusoids and attempt to recover those rotation periods.
    Saves an array of injected periods (days), recovered periods (days),
    rmag, injected amplitudes (ppm) and noise (ppm).
    """

    pers = np.exp(np.random.uniform(.1, 300, N))
    amps = np.random.uniform(1e2, 1e5, N)
    rmags = np.random.uniform(16, 28, N)
    noises_mag = np.array([LSST_sig(mag) for mag in rmags])
    noises_ppm = (1 - 10**(-noises_mag/2.5)) * 1e6

    # Simulate light curves
    print("Simulating light curves...")
    path = "simulations/sin"  # where to save the lcs
    [simulate_sin(i, b, pers[i], amps[i], path, noises_ppm[i]) for i in
     range(len(pers))]

    print("Saving results")
    data = np.vstack((pers, amps, rmags, noises_ppm))
    np.savetxt("sin_parameters_{0}.txt".format(fname), data.T)
    return pers, amps, rmags, noises_ppm

if __name__ == "__main__":
    fname = "sin_b{0}".format(sys.argv[1])

#     # Run simulations
    N = 30
    pers, amps, rmags, noises_ppm = inject("{0}".format(fname), N)

    # # recover periods
    # pers, amps, teffs, rmags, noises_ppm = \
    #         np.genfromtxt("parameters_{0}.txt".format(fname)).T
    # N = len(pers)
    # years = [1, 5, 10]
    # for year in years:
    #     periods = pgram(N, year, fname)
    #     data = np.vstack((pers, periods, np.log(amps), teffs, rmags, amps,
    #                       noises_ppm))
    #     np.savetxt("{0}yr_results{1}.txt".format(year, fname), data.T)

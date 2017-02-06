from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import mklc
from LSSToy import generate_visits
import os

plotpar = {'axes.labelsize': 20,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'text.usetex': True}
plt.rcParams.update(plotpar)

def simulate_LSST(x, depth, filt, id, p, a, path, noise, tmin=3, tmax=30,
                  dur=10, plot=False):
    """ Photometry with precision of 10 ppm (?).
    Uneven time sampling that ranges (uniformily) from 3 to 30 days (?).
    Lasting 10 years (?).
    id: id of the star.
    p: rotation period in seconds
    a: amplitude in ppm
    path: path to save files
    tmin, tmax: min and max intervals between observations in days.
    dur: duration in years.
    noise: noise level (ppm). Default is 10 ppm.
    """
    print(id)
    id = str(int(id)).zfill(4)

    sin2incl = np.random.uniform(np.sin(0)**2, np.sin(np.pi/2)**2)
    incl = np.arcsin(sin2incl**.5)
    tau = np.exp(np.random.uniform(np.log(p), np.log(10*p)))
    res0, res1 = mklc.mklc(x, incl=incl, tau=tau, p=p)
    # res0, res1 = mklc.mklc(x, p=p)
    nspot, ff, amp_err = res0
    time, area_tot, dF_tot, dF_tot0 = res1
    noise_free_y = dF_tot0 / np.median(dF_tot0) - 1
    y = dF_tot0 / np.median(dF_tot0) - 1 + noise*1e-6 * \
        np.random.randn(len(x))
    yerr = np.ones(len(y)) * noise * 1e-6

    data = np.vstack((x, y, yerr))
    np.savetxt("{0}/{1}.txt".format(path, id), data.T)
    # truths = np.array([p, a])
    # truths = np.vstack((np.array([p]), np.array([a])))
    # np.savetxt(f, truths.T)
    f = open("{0}/all_truths.txt".format(path), "a")
    f.write("{0} {1}\n".format(p, a))
    f.close()

    if plot:
        u = filt == "u"
        g = filt == "g"
        r = filt == "r"
        i = filt == "i"
        z = filt == "z"
        y = filt == "y"
        print("plotting light curve")
        plt.clf()
        plt.errorbar(x[u]/365.25, y[u], yerr=yerr[u], fmt=".", capsize=0,
                     ecolor=".5", color="b", label="u")
        plt.errorbar(x[g]/365.25, y[g], yerr=yerr[g], fmt=".", capsize=0,
                     ecolor=".5", color="g", label="g")
        plt.errorbar(x[r]/365.25, y[r], yerr=yerr[r], fmt=".", capsize=0,
                     ecolor=".5", color="r", label="r")
        plt.errorbar(x[i]/365.25, y[i], yerr=yerr[i], fmt=".", capsize=0,
                     ecolor=".5", color="m", label="i")
        plt.errorbar(x[z]/365.25, y[z], yerr=yerr[z], fmt=".", capsize=0,
                     ecolor=".5", color="y", label="z")
        plt.errorbar(x[y]/365.25, y[y], yerr=yerr[y], fmt=".", capsize=0,
                     ecolor=".5", color="k", label="y")
        plt.xlabel("$\mathrm{Time~(years)}$")
        plt.ylabel("$\mathrm{Normalised~Flux}$")
        plt.xlim(min(x/365.25), max(x/365.25))
        plt.subplots_adjust(left=.2, bottom=.12)
        plt.savefig(os.path.join(path, "{}".format(id)))

if __name__ == "__main__":
    path = "simulations"  # where to save

    # Arrays of random (log-normal) periods and (uniform) amplitudes.
    N = 10
    ps = np.exp(np.random.uniform(np.log(2), np.log(100), N))
    amps = np.random.uniform(10, 300, N)  # ppm
    [simulate_LSST(i, ps[i], amps[i], path) for i in range(N)]

    # save the true values
    ids = np.arange(N)
    data = np.vstack((ids, ps, amps))
    np.savetxt("{0}/truth.txt".format(path), data.T)

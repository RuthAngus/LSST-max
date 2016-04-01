from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import mklc
from LSSToy import generate_visits
import os

def simulate_LSST(id, p, a, path, noise, tmin=3, tmax=30, dur=10):
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
    id = str(int(id)).zfill(4)

    # only make new simulation if one doesn't already exist
    if os.path.exists("simulations/{0}".format(id)):
        return

    # The time array
#     x = np.cumsum(np.random.uniform(tmin, tmax, 1000))
#     x = x[x < dur * 365.25]
    x = generate_visits()
    x += -x[0]

    np.random.seed(1234)
    res0, res1 = mklc.mklc(x, p=p)
    nspot, ff, amp_err = res0
    time, area_tot, dF_tot, dF_tot0 = res1
    y = dF_tot0 / np.median(dF_tot0) - 1 + noise*1e-6 * np.random.randn(len(x))
    yerr = np.ones_like(y) * noise * 1e-6

    data = np.vstack((x, y, yerr))
    np.savetxt("{0}/{1}.txt".format(path, id), data.T)

    plt.clf()
    plt.errorbar(x/365.25, y, yerr=yerr, fmt="k.", capsize=0)
    plt.xlabel("Time (years)")
    plt.ylabel("Normalised flux")
    plt.title("period = {0:.2f} days, amp = {1:.2f} ppm".format(p, a))
    plt.subplots_adjust(left=.2)
    plt.savefig("simulations/{0}".format(id))

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

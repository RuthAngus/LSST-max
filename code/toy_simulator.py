from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as spi
import mklc
# from simlc import spots

def suz():
    # Load one of Suzanne's simulated kepler light curves
    spath = "../../GProtation/code/simulations/kepler_diffrot_full/"
    data = np.genfromtxt("{0}/par/final_table.txt".format(spath),
                         skip_header=1).T
    m = data[13] == 0  # just the stars without diffrot
    ids = data[0][m][:2]
    xs, ys = np.genfromtxt("{0}/noise_free/lightcurve_{1}.txt".format(spath,
                                                                      id)).T
    f = spi.interp1d(xs, ys, kind="cubic")
    y = f(x)  # interpolate onto LSST cadence
    y += noise/1e6 * np.random.randn(len(x))
    yerr = np.ones_like(y) * noise/1e6

def simulate_LSST(id, p, a, path, tmin=3, tmax=30, dur=10, noise=10.):
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

#     xs = np.arange(0, 10*365.25)
#     np.random.seed(1234)
#     res0, res1 = mklc.mklc(xs, p=p)
#     nspot, ff, amp_err = res0
#     time, area_tot, dF_tot, dF_tot0 = res1
#     ys = dF_tot0 / np.median(dF_tot0) - 1

    # The time array
    x = np.cumsum(np.random.uniform(tmin, tmax, 1000))
    x = x[x < dur * 365.25]  # cut off at 10 yrs

    np.random.seed(1234)
    res0, res1 = mklc.mklc(x, p=p)
    nspot, ff, amp_err = res0
    time, area_tot, dF_tot, dF_tot0 = res1
    y = dF_tot0 / np.median(dF_tot0) - 1 + noise*1e-6 * np.random.randn(len(x))
    yerr = np.ones_like(y) * noise/1e6

#     spots = spots()
#     area, ome, beta, dF = spots.calc(

    # The flux array (for now just use a sinusoid)
#     y = (float(a)/1e6) * np.sin(2*np.pi*(1./p)*x) \
#             + noise/1e6 * np.random.randn(len(x))
#     yerr = np.ones_like(y) * noise/1e6
#     xs = np.linspace(min(x), max(x), 1000)
#     ys = (float(a)/1e6) * np.sin(2*np.pi*(1./p)*xs)

    data = np.vstack((x, y, yerr))
    np.savetxt("{0}/{1}.txt".format(path, id), data.T)

    plt.clf()
#     plt.plot(xs/365.25, ys, "r")
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
    ps = np.exp(np.random.uniform(np.log(10), np.log(300), N))
    amps = np.random.uniform(10, 300, N)  # ppm
    [simulate_LSST(i, ps[i], amps[i], path) for i in range(N)]

    # save the true values
    ids = np.arange(N)
    data = np.vstack((ids, ps, amps))
    np.savetxt("{0}/truth.txt".format(path), data.T)

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

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

    # The time array
    x = np.cumsum(np.random.uniform(tmin, tmax, 1000))
    x = x[x < dur * 365.25]  # cut off at 10 yrs

    # The flux array (for now just use a sinusoid)
    y = (float(a)/1e6) * np.sin(2*np.pi*(1./p)*x) \
            + noise/1e6 * np.random.randn(len(x))
    yerr = np.ones_like(y) * noise/1e6
    xs = np.linspace(min(x), max(x), 1000)
    ys = (float(a)/1e6) * np.sin(2*np.pi*(1./p)*xs)

    data = np.vstack((x, y, yerr))
    np.savetxt("{0}/{1}.txt".format(path, id), data.T)

    plt.clf()
    plt.plot(xs/365.25, ys, "r")
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

import numpy as np
import matplotlib.pyplot as plt
from gatspy.periodic import LombScargle
import scipy.signal as sps

def window(b):

    time, depth = np.genfromtxt("data/l45b{0}_cadence.txt".format(b)).T
    ps = np.linspace(.1, 100, 1000)
    fs = 1./ps
    ys = np.ones_like(time)
    normval = time.shape[0]
    pgram = sps.lombscargle(time, ys, 2*np.pi*fs)
    return ps, pgram, normval, time, depth

if __name__ == "__main__":
    bs = [-10.0, -20.0, -40.0, -60.0, -80.0]
    plt.clf()
    for b in bs:
        ps, pgram, normval, t, d = window(b)
        plt.plot(ps, np.sqrt(4*(pgram/normval)), label=b)
    plt.legend()
    plt.savefig("window.pdf".format(b))

    plt.clf()
    for b in bs:
        ps, pgram, normval, t, d = window(b)
        plt.plot(t, d, ".", label=b)
    plt.legend()
    plt.savefig("cadences_for_all_bs.pdf".format(b))

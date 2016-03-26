import numpy as np
import matplotlib.pyplot as plt
import barnes as bn


def gr_to_bv(g, r):
    """
    From https://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
    Returns the B-V value and the RMS error.
    """
    return .62*(g - r) + .15, .07


def random_stars(fname, N):
    """
    Randomly draw stellar properties from the output of a trilegal field.
    fname: str, the name of the trilegal output file, e.g.
    "output574523944248.dat"
    N: int, the number of stars you want to draw.
    Returns an array of ages and B-V colours
    """
    # load data
    Gc, logAge, m_h, m_ini, logL, logTeff, logg, m_M0, Av, m2_m1, mbol, u, g,\
        r, i, z, Mact = np.genfromtxt(fname).T

    # convert to B-V and get rid of hot stars
    bv, bverr = gr_to_bv(g, r)
    # m = bv > .4  # no hot stars
    # logAge, g, r = logAge[m], g[m], r[m]

    # randomly select stars.
    stars = np.random.choice(np.arange(len(logAge)), N)
    logAges, gs, rs = logAge[stars], g[stars], r[stars]
    bvs, bverrs = gr_to_bv(gs, rs)
    return logAges, bvs


if __name__ == "__main__":
    fname = "output574523944248.dat"
    nstars = 1000
    logAges, bvs = random_stars(fname, nstars)

    age = 10**logAges * 1e-6
    ps = bn.period(age, bvs)

    plt.clf()
    plt.hist(ps)
    plt.savefig("period_hist")

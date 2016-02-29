from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import h5py
from plotstuff import params, colours
reb = params()
cols = colours()
from gatspy.periodic import LombScargle
import sys
import multiprocessing as mp
from multiprocessing import Pool
from GProtation import make_plot, lnprob
import emcee
import time
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel

def periodograms(id, x, y, yerr, path, plot=False, savepgram=False):
    """
    takes id of the star, returns an array of period measurements and saves the
    results.
    id: star id.
    x, y, yerr: time, flux and error arrays.
    path: path where you want to save the output.
    """
    ps = np.linspace(5, 300, 1000)
    model = LombScargle().fit(x, y, yerr)
    pgram = model.periodogram(ps)

    # find peaks
    peaks = np.array([i for i in range(1, len(ps)-1) if pgram[i-1] <
                     pgram[i] and pgram[i+1] < pgram[i]])
    if len(peaks):
        period = ps[pgram==max(pgram[peaks])][0]
    else: period = 0

    if plot:
        plt.clf()
        plt.plot(ps, pgram)
        plt.axvline(period, color="r")
        plt.savefig("{0}/{1}_pgram".format(path, str(int(id)).zfill(4)))

    if savepgram:
        np.savetxt("{0}/{1}_pgram.txt".format(path, str(int(id)).zfill(4)),
                   np.transpose((ps, pgram)))

    np.savetxt("{0}/{1}_pgram_result.txt".format(path, str(int(id)).zfill(4)),
               np.ones(2).T*period)
    return period

def recover_injections(id, x, y, yerr, path, burnin, run, nwalkers=32,
                       plot=True):
    """
    Take x, y, yerr, calculate ACF period for initialisation and do MCMC.
    npts: number of points per period.
    id: star id.
    x, y, yerr: time, flux and error arrays.
    path: path where you want to save the output.
    burnin: the number of burnin steps.
    run: the number of steps to run for.
    nwalkers: the number of walkers.
    plot: if True then plots of posteriors and chains will be made.
    """

    # initialise with pgram
    try:
        p_init = np.genfromtxt("{0}/{1}_pgramresult.txt".format(path, id))
    except:
        p_init = periodograms(id, x, y, yerr, path, plot=True)

    if p_init < .5:  # small periods raise an error with george.
            p_init = 1.

    # Set limits on prior
    plims = np.log([p_init - .4 * p_init, p_init + .4 * p_init])

    print("Initial period and limits:", p_init, np.exp(plims))

    # assign theta_init
    theta_init = [np.exp(-5), np.exp(7), np.exp(.6), np.exp(-16), p_init]
    theta_init = np.log(theta_init)
    print("\n", "log(theta_init) = ", theta_init)
    print("theta_init = ", np.exp(theta_init), "\n")

    # set up MCMC
    ndim, nwalkers = len(theta_init), nwalkers
    p0 = [theta_init+1e-4*np.random.rand(ndim) for i in range(nwalkers)]
    args = (x, y, yerr, plims)

    # time the lhf call
    start = time.time()
    print("lnprob = ", lnprob(theta_init, x, y, yerr, plims))
    end = time.time()
    tm = end - start
    print("1 lhf call takes ", tm, "seconds")
    print("burn in will take", tm * nwalkers * burnin, "s")
    print("run will take", tm * nwalkers * run, "s")
    print("total = ", (tm*nwalkers*run + tm*nwalkers*burnin)/60, \
          "mins,", (tm*nwalkers*run + tm*nwalkers*burnin)/3600, "hours")

    # run MCMC
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=args)
    print("burning in...")
    start = time.time()
    p0, lp, state = sampler.run_mcmc(p0, burnin)
    sampler.reset()
    print("production run...")
    p0, lp, state = sampler.run_mcmc(p0, run)
    end = time.time()
    print("actual time = ", (end - start)/60, "mins")

    # save samples
    f = h5py.File("%s/%s_samples.h5" % (path, id), "w")
    data = f.create_dataset("samples", np.shape(sampler.chain))
    data[:, :] = np.array(sampler.chain)
    f.close()

    # make various plots
    if plot:
        with h5py.File("%s/%s_samples.h5" % (path, id), "r") as f:
            samples = f["samples"][...]
        mcmc_result = make_plot(samples, x, y, yerr, id, path, traces=True,
                                tri=True, prediction=True)

def acf_pgram_GP_LSST(id):
    """
    Run acf, pgram and MCMC recovery on Suzanne's simulations
    """
    id = str(int(id)).zfill(4)
    path = "results"  # where to save results
    x, y, yerr = np.genfromtxt("simulations/{0}.txt".format(id)).T

    periodograms(id, x, y, yerr, path, plot=True)  # pgram
    burnin, run = 50, 1000
    recover_injections(id, x, y, yerr, path, burnin, run, nwalkers=32,
                       plot=True)

if __name__ == "__main__":

    acf_pgram_GP_LSST(0)

#     ids = range(10)
#     pool = Pool()
#     pool.map(acf_pgram_GP_LSST, ids)

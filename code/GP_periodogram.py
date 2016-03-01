from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from GProtation import make_plot, lnprob, neglnlike
import emcee
import time
import george
from george.kernels import ExpSquaredKernel, ExpSine2Kernel
import scipy.optimize as spo

def GP_periodogram(x, y, yerr, p_init, plims, N):
	"""
	This function takes a light curves and attempts to produce a GP periodogram.
	It returns the value of the highest peak.
	The kernel hyperparameters are optimised over a grid of periods.
	This is also a "profile likelihood".
	x, y, yerr: the light curve.
	p_init: the initial guess for the period.
	plims: the (log) boundaries for the grid.
	N: the number of grid points.
	"""

	# create the grid
	periods = np.linspace(np.exp(plims[0], np.exp(plims[1], 10)

	# initial hyperparameters

if __name__ == "__main__":

	# fake data
	x = np.arange(0, 10, 100)
	p = 2
	err = .1
	y = np.sin(2*np.pi*(1./p)*x) + np.random.randn(100)*err
	yerr = np.ones_like(y) * err
	p_init, plims = 2, np.log(.1, 5)

	GP_periodogram(x, y, yerr, p_init, plims, 10)

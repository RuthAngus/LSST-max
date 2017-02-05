# Calculating the star spot variations for the Wfirst bandpass.

import os

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as si

from astropy.constants import c, h, k_B


def planck(lda, T):
    """
    Calculate the spectral radiance of a star at a given temperature.
    params:
    ------
    lda: (float)
        Wavelength to query.
    T: (float)
        Temperature of star.
    """
    exponent = h*c/(lda*k_B*T)
    spectral_radiance = (2*h*c**2/lda**5) * 1./(np.exp(exponent.value) - 1.)
    return spectral_radiance.value


def relative_fluxes(lower, upper, starteff, spotteff, fraction):
    """
    Total relative flux in the Wfirst main, IR bandpass for a given percentage
    spot coverage.
    params:
    ------
    lower: (float)
        lower wavelength of bandpass.
    lower: (float)
        lower wavelength of bandpass.
    starteff: (float)
        Temperature of star.
    spotteff: (float)
        Temperature of spot.
    fraction: (float)
        The fraction of the visible stellar surface taken up by a spot.
    Returns:
    -------
    frac_variation: (float)
        fractional flux variation due to spots for a given star teff, spot
        teff and bandpass.

    """
    # Total unspotted star flux.
    flux = si.quad(planck, lower, upper, args=(starteff))[0]

    # Flux of star-sized spot.
    spotflux = si.quad(planck, lower, upper, args=(spotteff))[0]

    # Flux of spotted star.
    totflux = flux - flux*fraction + spotflux*fraction

    # Difference between spotted star flux and non-spotted star flux
    variation = flux - totflux
    frac_variation = variation/flux
    return frac_variation


if __name__ == "__main__":
    wavelengths = np.linspace(1e-7, 5e-6, 1000)
    teff = 5000
    teff_diff = 2000
    spot_teff = teff - teff_diff
    star = planck(wavelengths, 5000)
    starspot = planck(wavelengths, 2000)

    plt.clf()
    plt.plot(wavelengths, star)
    plt.plot(wavelengths, starspot)
    plt.savefig("test")

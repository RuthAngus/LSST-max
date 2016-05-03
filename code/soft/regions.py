import numpy as np
import matplotlib.pyplot as plt
import time

def regions(seed=0, randspots, activityrate=1, cyclelength=1, cycleoverlap=0,
            maxlat=70, minlat=0, tsim=1000, tstart=0, dir="."):
    """
    inputs
    activityrate - number of bipoles (1= solar)
    cyclelength - length of cycle in years
    cycleoverlap - cycleoverlap time in years
    tsim - length of simulation in days
    tstart - first day to start outputting bipoles
    minlat - minimum latitude of spot emergence
    maxlat - maximum latitude of spot emergence
    randspots - set with /randspots. Use this for no cycle

    This program simulates the solar cycle. It produces a list
    of active regions with the following parameters:

    nday = day of emergence
    thpos= theta of positive pole (radians)
    phpos= phi   of positive pole (radians)
    thneg= theta of negative pole (radians)
    phneg= phi   of negative pole (radians)
    width= width of each pole (radians)
    bmax = maximum flux density (Gauss)

    According to Schrijver and Harvey (1994), the number of active regions
    emerging with areas in the range [A,A+dA] in a time dt is given by

    n(A,t) dA dt = a(t) A^(-2) dA dt ,

    where A is the "initial" area of a bipole in square degrees, and t is
    the time in days; a(t) varies from 1.23 at cycle minimum to 10 at cycle
    maximum.

    The bipole area is the area within the 25-Gauss contour in the
    "initial" state, i.e. time of maximum development of the active region.
    The assumed peak flux density in the initial sate is 1100 G, and
    width = 0.2*bsiz (see disp_region). The parameters written onto the
    file are corrected for further diffusion and correspond to the time
    when width = 4 deg, the smallest width that can be resolved with lmax=63.

    In our simulation we use a lower value of a(t) to account for "correlated"
    regions.
    """
    nbin = 5                              # number of area bins
    delt = 0.5                            # delta ln(A)
    amax = 100.                           # orig. area of largest bipoles (deg^2)
    dcon = exp(0.5*delt)-exp(-0.5*delt)   # contant from integ. over bin

    print('Creating regions with the following parameters:')
    print('Activity rate: ', activityrate, ' x Solar rate.'
    print('Cycle length: ', cyclelength, ' years.')
    print('Cycle overlap: ', cycleoverlap, ' years.')
    print('Max spot lat: ', maxlat, ' degrees')
    print('Min spot lat: ', minlat, ' degrees')
    print('Simulation time: ', tsim, ' days')
    print('Simulation start:', tstart, ' days')
    deviation = 5   # rmsd deviation from butterfly pattern
    atm = np.zeros(200) + 10.*activityrate
    ncycle  = np.zeros(200) + cyclelength
    nclen   = np.zeros(200) + cyclelength + cycleoverlap
    latrmsd = np.zeros(200) + deviation

    # a(t) at cycle maximum (deg^2/day)
    # cycle period (days)
    # cycle duration (days)

    ncycle, nclen = ncycle*365., nclen*365.
    fact = np.exp(delt*np.ones(nbin))  # array of area reduction factors
    ftot = sum(fact)  # sum of reduction factors
    bsiz = np.sqrt(amax/fact)  # array of bipole separations (deg)
    tau1 = 5.  # first and last times (in days) for
    tau2 = 15.  # emergence of "correlated" regions
    prob = 0.001  # total probability for "correlation"
    nlon = 36  # number of longitude bins
    nlat = 16  # number of latitude bins
    nday1 = 0  # first day to be simulated
    ndays = tsim  # number of days to be simulated
    dt = 1


def add_region(nday, ic, lo, lat, k, bsiz1, seed, phase):
    w_org = .4 * bsiz1  # original width (degrees), at birth
    width = 4.  # final width (degrees), at death
    bmax = 250. * (w_org / width)**2  # final peak flux density (G)
    bsizr = np.pi * bsiz1 / 180.  # pole separation in radians
    width = np.pi * width / 180.  # final width in radians

    np.seed(seed)
    rand_array = np.random.randn(100)  # random number less than 1.6
    x = rand_array[rand_array < 1.6][0]
    y = rand_array[rand_array < 1.8][0]
    np.seed(seed)
    z = np.random.uniform(1)

if __name__ == "__main__":


    # # Initialize random number generator:
    # if seed == -1:
    #       seed = time.time()

    # # Initialize time since last emergence of a large region, as function
    # # of longitude, latitude and hemisphere:

    # tau = range(nlon, nlat, 2) + tau2
    # dlon = 360. / nlon
    # dlat = maxlat / nlat

    # ncnt = 0
    # # Loop over time (in days):

    # ncur = 0
    # cycle_days = ncycle[0]
    # start_day = 0

    # for nday in range(ndays):

    #     # Compute index of most recently started cycle:
    #     ncur_test = nday / cycle_days
    #     ncur_now = nday / cycle_days
    #     ncur_prev = (nday - 1) / cycle_days

    #     if int(ncur_now) != int(ncur_prev):
    #         ncur += ncur
    #         cycle_days += ncycle(ncur)

    # # Initialize rate of emergence for largest regions, and add 1 day
    # # to time of last emergence:

    # tau += 1
    # rc0 = np.arange(nlon, nlat, 2)
    # index = (tau > tau1) * (tau l<= tau2)
    #     if index[0] > -1:
    #         rc0[index] = prob / (tau2 - tau1)

    # # Loop over current and previous cycle:
    # for icycle in range(2):
    #     nc = ncur - icycle  # index of cycle
    #     if ncur == 0:
    #         nc1 = 0
    #         start_day = nc*ncycle[0]
    #     else:
    #         nc1 = nc
    #         if ncur == 1:
    #             if icycle == 0:
    #                 start_day = int(sum(ncycle[:nc-1]))
    #                 if icycle == 1:
    #                     start_day = 0
    #         else:
    #             start_day = int(sum(ncycle[:nc-1]))

    # np.savetxt("{0}/regions.txt".format(dir))


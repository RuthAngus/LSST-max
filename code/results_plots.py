import numpy as np
import matplotlib.pyplot as plt
from plotstuff import colours
cols = colours()

plotpar = {'axes.labelsize': 16,
           'xtick.labelsize': 16,
           'ytick.labelsize': 16,
           'text.usetex': True}
plt.rcParams.update(plotpar)


def make_plot(data, b):
    """
    Scatter plots of recovered period vs injected periods, coloured by
    r-magnitude, teff and amplitude.
    """

    pers, periods, log_amps, teffs, rmags, amps, noises_ppm = data

    xs = np.linspace(0, max(pers))

    plt.clf()
    plt.plot(xs, xs, ":", color="k")
    plt.scatter(pers, periods, c=teffs, vmin=3000, vmax=8000,
                edgecolor="face", cmap="BuPu", s=8)
    plt.colorbar(label="$\mathrm{T}_{\mathrm{eff}}~\mathrm{(K)}$")
    plt.ylabel("$\mathrm{Measured~Period~(Days)}$")
    plt.xlabel("$\mathrm{Injected~Period~(Days)}$")
    plt.xlim(0, max(pers))
    plt.ylim(0, max(periods))
    plt.savefig("pvp_T_{0}.pdf".format(b))

    plt.clf()
    plt.plot(xs, xs, ":", color="k")
    plt.scatter(pers, periods, c=rmags,
                edgecolor="face", cmap="GnBu_r", s=8)
    plt.colorbar(label="$\mathrm{r-band~magnitude}$")
    plt.ylabel("$\mathrm{Measured~Period~(Days)}$")
    plt.xlabel("$\mathrm{Injected~Period~(Days)}$")
    plt.xlim(0, max(pers))
    plt.ylim(0, max(periods))
    plt.savefig("pvp_r_{0}.pdf".format(b))

    plt.clf()
    plt.plot(xs, xs, ":", color="k")
    plt.scatter(pers, periods, c=log_amps,
                edgecolor="face", cmap="PuRd", s=8)
    plt.colorbar(label="$\log\mathrm{(Amplitude)~(ppt)}$")
    plt.ylabel("$\mathrm{Measured~Period~(Days)}$")
    plt.xlabel("$\mathrm{Injected~Period~(Days)}$")
    plt.xlim(0, max(pers))
    plt.ylim(0, max(periods))
    plt.savefig("pvp_a_{0}.pdf".format(b))


def recovered(data, f):
    """
    Take all data and return just the successfully recovered stuff
    data: 2d array of pers, periods, log_amps, teffs, rmags, amps, noises_ppm
    f: tolerance, e.g. 0.1.
    Returns 2d array of just successfully recovered stars.
    """
    pers, periods, log_amps, teffs, rmags, amps, noises_ppm = data
    m = (periods < pers + f*pers) * (pers - f*pers < periods)
    return np.vstack((pers[m], periods[m], log_amps[m], teffs[m], rmags[m],
                      amps[m], noises_ppm[m]))


def find_fraction(X, X_r, bins):
    # Fraction recovered as a function of X
    if bins is None:
        true_hist, bins = np.histogram(X)
        measured_hist, _ = np.histogram(X_r, bins)
    else:
        true_hist, _ = np.histogram(X, bins)
        measured_hist, _ = np.histogram(X_r, bins)
    th = np.array([float(i) for i in true_hist])
    mh = np.array([float(i) for i in measured_hist])
    return mh/th*100, bins  # percent


def percents(X, X_r, teffs, teffs_r):
    # total --- all stars
    percent, bins = find_fraction(X, X_r, None)

    # now for different temperatures
    mf = lambda Xmin, Xmax, t: (Xmin < t) * (t < Xmax)
    Gm, Km, Mm, Fm = mf(Gmin, Gmax, teffs), mf(Kmin, Kmax, teffs), \
        mf(Mmin, Mmax, teffs), mf(Fmin, Fmax, teffs)
    Gmr, Kmr, Mmr, Fmr = mf(Gmin, Gmax, teffs_r), mf(Kmin, Kmax, teffs_r), \
        mf(Mmin, Mmax, teffs_r), mf(Fmin, Fmax, teffs_r)
    Gpercent, Gbins = find_fraction(X[Gm], X_r[Gmr], bins)
    Kpercent, Kbins = find_fraction(X[Km], X_r[Kmr], bins)
    Mpercent, Mbins = find_fraction(X[Mm], X_r[Mmr], bins)
    Fpercent, Fbins = find_fraction(X[Fm], X_r[Fmr], bins)
    return Gpercent, Kpercent, Mpercent, Fpercent, bins


def summary_plot(data, b, f):
    """
    Make plots of completeness vs rotation period.
    Completeness vs amplitude.
    Completeness vs r-mag.
    Completeness vs teff.
    """

    # load data
    pers, periods, log_amps, teffs, rmags, amps, noises_ppm = data

    # find the arrays of successful recoveries
    pers_r, periods_r, log_amps_r, teffs_r, rmags_r, amps_r, noises_ppm_r = \
        recovered(data, f)
    print("\n", len(pers), "injected", len(pers_r),
          "recovered, {0:.2f}%".format(float(len(pers_r)) /
                                       float(len(pers))*100), "\n")

    # Calculate the percentage recovery as a function of period for the
    # different spectral types
    Gpercent, Kpercent, Mpercent, Fpercent, bins = percents(pers, pers_r,
                                                            teffs, teffs_r)

    Gfrac, Kfrac, Mfrac = 0.2426, 0.5429, 0.8339
    # Gfrac, Kfrac, Mfrac = 1, 1, 1
    # make the plot
    plt.clf()
    # plt.step(bins[:-1], Fpercent, lw=2, color="MediumPurple",
    #          label="$\mathrm{F~dwarfs}$")
    plt.step(bins[:-1], Gpercent*Gfrac, lw=2, color="CornflowerBlue",
             label="$\mathrm{G~dwarfs}$")
    plt.step(bins[:-1], Kpercent*Kfrac, lw=2, color="LimeGreen",
             label="$\mathrm{K~dwarfs}$")
    plt.step(bins[:-1], Mpercent*Mfrac, lw=2, color="DarkOrange",
             label="$\mathrm{M~dwarfs}$")
    plt.xlim(0, bins[-2])
    plt.xlabel("$\mathrm{Injected~Rotation~Period~(Days)}$")
    plt.ylabel("$\mathrm{Percentage~Successfully~Recovered}$")
    plt.ylim(0, 100)
    plt.legend(loc="best")
    plt.savefig("recovered_hist_{0}.pdf".format(b))


def trilegal(data, b):
    """
    Make histograms of trilegal outputs.
    Histograms of Rotation period for different masses
    """

    # load data
    pers, periods, log_amps, teffs, rmags, amps, noises_ppm = data

    mf = lambda Xmin, Xmax, t: (Xmin < t) * (t < Xmax)

    m = (teffs < 7000) * (2000 < teffs)

    plt.clf()
    plt.hist(teffs[m], color="w", histtype="stepfilled")
    plt.xlabel("$\mathrm{Effective~Temperature~(K)}$")
    plt.ylabel("$\mathrm{Number~out~of~20,0000}$")
    plt.savefig("trilegal_teff_hist{0}.pdf".format(b))

    plt.clf()
    m = (teffs < 7000) * (2000 < teffs) * pers > 0
    bins = plt.hist(pers[m], 10, color="w", histtype="stepfilled")[1]
    phist = np.histogram(pers[m], bins=bins)[0]
    plt.hist(pers[m], bins=bins, color="w", histtype="stepfilled")
    rec = recovered(data, .1)
    l = (rec[3] < 7000) * (2000 < rec[3]) * rec[0] > 0
    plt.hist(rec[0][l], bins=bins, color="w", edgecolor="CornFlowerBlue",
             histtype="stepfilled", lw=2)
    plt.xlabel("$\mathrm{Injected~Rotation~Period~(Days)}$")
    plt.ylabel("$\mathrm{Number~out~of~20,0000}$")
    plt.savefig("trilegal_period_hist{0}.pdf".format(b))

    # make fraction recovered histogram
    nbins = 10
    m = (teffs < 7000) * (2000 < teffs) * pers > 0
    rec = recovered(data, .1)
    l = (rec[3] < 7000) * (2000 < rec[3]) * rec[0] > 0

    # make histograms
    inj, ibins = np.histogram(pers[m], nbins)
    reco, rbins = np.histogram(rec[0][l], bins=ibins)
    frac = np.zeros(len(ibins))
    frac[:-1] = reco/inj
    kepler_frac = 0.328
    plt.clf()
    plt.step(ibins, frac*kepler_frac, color="CornFlowerBlue", lw=2)
    # plt.xlim(0, 70)
    # plt.ylim(0, 1)
    plt.xlabel("$\mathrm{Injected~Rotation~Period~(Days)}$")
    plt.ylabel("$\mathrm{Fraction~Recovered}$")

    # load Amy's data
    data = np.genfromtxt("Table_1_Periodic.txt", skip_header=1,
                         delimiter=",").T
    m = (2000 < data[1]) * (data[1] < 7000)
    amy_p = data[4][m]
    amy, abins = np.histogram(amy_p, bins=ibins)
    amy_tot = 103693. / nbins
    # plt.step(ibins[:-1], amy*(phist/float(max(phist))), "k")
    # plt.step(ibins[:-1], amy/amy_tot/phist, "k")
    phist = np.array([float(i) for i in phist])
    print(amy, phist/float(max(phist)))
    # plt.step(ibins[:-1], amy*phist, "k")
    # plt.step(bins[:-1], phist/float(max(phist)), "r")
    # plt.step(ibins[:-1], amy/amy_tot * phist/float(max(phist)), "k")
    plt.ylim(0, .5)
    plt.savefig("trilegal_fraction_amy_hist{0}.pdf".format(b))


if __name__ == "__main__":

    Fmin, Fmax = 6000, 7000
    Gmin, Gmax = 5200, 6000
    Kmin, Kmax = 3700, 5200
    Mmin, Mmax = 2400, 3700

    yr = 10
    b = -10

    # load data file
    d = np.genfromtxt("results/{0}yr_resultsl45b{1}.txt".format(yr, b)).T
    m = d[0] > 0  # remove negative injected periods
    data = np.vstack((d[0][m], d[1][m], d[2][m], d[3][m], d[4][m], d[5][m],
                      d[6][m]))

    # make the plots
    # summary_plot(data, b, .1)
    # make_plot(data, b)
    trilegal(data, b)

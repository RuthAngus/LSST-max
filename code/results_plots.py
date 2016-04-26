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
    """
    pers, periods, log_amps, teffs, rmags, amps, noises_ppm = data
    m = (periods < pers + f*pers) * (pers - f*pers < periods)
    return np.vstack((pers[m], periods[m], log_amps[m], teffs[m], rmags[m],
                      amps[m], noises_ppm[m]))

def find_fraction(pers, pers_r):
    # Fraction recovered as a function of P
    true_hist, bins = np.histogram(pers)
    measured_hist, _ = np.histogram(pers_r, bins)
    th = np.array([float(i) for i in true_hist])
    mh = np.array([float(i) for i in measured_hist])
    return mh/th*100, bins  # percent

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
          "recovered, {0:.2f}%".format(float(len(pers_r))/
                                       float(len(pers))*100), "\n")

    # Calculate the percentage recovery as a function of period for the
    # different spectral types
    Gmin, Gmax = 5200, 6000
    Kmin, Kmax = 3700, 5200
    Mmin, Mmax = 2400, 3700
    mf = lambda Xmin, Xmax, t: (Xmin < t) * (t < Xmax)
    Gm, Km, Mm = mf(Gmin, Gmax, teffs), mf(Kmin, Kmax, teffs), \
                 mf(Mmin, Mmax, teffs)
    Gmr, Kmr, Mmr = mf(Gmin, Gmax, teffs_r), mf(Kmin, Kmax, teffs_r), \
                 mf(Mmin, Mmax, teffs_r)
    Gpercent, Gbins = find_fraction(pers[Gm], pers_r[Gmr])
    Kpercent, Kbins = find_fraction(pers[Km], pers_r[Kmr])
    Mpercent, Mbins = find_fraction(pers[Mm], pers_r[Mmr])

    # make the plot
    plt.clf()
    plt.step(Gbins[:-1], Gpercent, lw=2, color=cols.blue)
    plt.step(Kbins[:-1], Kpercent, lw=2, color=cols.orange)
    plt.step(Mbins[:-1], Mpercent, lw=2, color=cols.pink)
    plt.xlim(0, Gbins[-2])
    plt.xlabel("$\mathrm{Injected~Rotation~Period~(Days)}$")
    plt.ylabel("$\mathrm{Percentage~Successfully~Recovered}$")
    plt.savefig("recovered_hist.pdf")

if __name__ == "__main__":
    yr = 10
    b = -10

    # load data file
    data = np.genfromtxt("{0}yr_resultsl45b{1}.txt".format(yr, b)).T
    m = data[0] > 0
    data2 = np.vstack((data[0][m], data[1][m], data[2][m], data[3][m],
                       data[4][m], data[5][m], data[6][m]))

    # make the plots
    summary_plot(data2, b, .1)
    # make_plot(data, b)

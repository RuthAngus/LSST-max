import numpy as np
import matplotlib.pyplot as plt
from plotstuff import colours
cols = colours()
plotpar = {'axes.labelsize': 20,
           'text.fontsize': 10,
           'legend.fontsize': 20,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20,
           'text.usetex': True}
plt.rcParams.update(plotpar)

fname = "output574523944248.dat"
fname = "output16533990464.dat"
data = np.genfromtxt("rotation_results{0}.txt".format(fname)).T

m = data[0] > 0
true = data[0][m]
recovered = data[1][m]

plt.clf()
xs = np.arange(100)
plt.plot(true, recovered, "k.", zorder=1)
plt.plot(xs, xs, color=cols.lightblue, lw=2, zorder=2)
# plt.plot(xs, xs, color=cols.blue, ls="--", zorder=2)
plt.fill_between(xs, xs - .1*xs, xs + .1*xs, color=cols.lightblue, alpha=.3,
                 zorder=0)
plt.fill_between(xs, xs - .2*xs, xs + .2*xs, color=cols.lightblue, alpha=.2,
                 zorder=0)
plt.xlim(0, 75)
plt.ylim(0, 100)
plt.xlabel("$\mathrm{True~Period~(days)}$")
plt.ylabel("$\mathrm{Measured~Period~(days)}$")
plt.subplots_adjust(bottom=.12)
plt.savefig("test.pdf")

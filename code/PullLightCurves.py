# coding: utf-8
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.db as db
import lsst.sims.maf.utils as utils
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles


def get_cadence(ra, dec, b, snrLimit, nPtsLimit, filters, outDir, opsimdb,
                resultsDb):

    # The pass metric just passes data straight through.
    metric = metrics.PassMetric(cols=['filter', 'fiveSigmaDepth', 'expMJD'])
    slicer = slicers.UserPointsSlicer(ra, dec, lonCol='ditheredRA',
                                      latCol='ditheredDec')
    sql = ''
    bundle = metricBundles.MetricBundle(metric, slicer, sql)
    bg = metricBundles.MetricBundleGroup({0: bundle}, opsimdb, outDir=outDir,
                                         resultsDb=resultsDb)

    bg.runAll()
    bundle.metricValues.data[0]['filter']

    print("Plotting...")
    colors = {'u': 'cyan', 'g': 'g', 'r': 'y', 'i': 'r', 'z': 'm', 'y': 'k'}
    dayZero = bundle.metricValues.data[0]['expMJD'].min()
    times = []
    depths = []
    plt.clf()
    for fname in filters:
        good = np.where(bundle.metricValues.data[0]['filter'] == fname)
        times.append(bundle.metricValues.data[0]['expMJD'][good] - dayZero)
        depths.append(bundle.metricValues.data[0]['fiveSigmaDepth'][good])

        plt.scatter(bundle.metricValues.data[0]['expMJD'][good] - dayZero,
                    bundle.metricValues.data[0]['fiveSigmaDepth'][good],
                    c=colors[fname], label=fname)

    plt.xlabel('Day')
    plt.ylabel('5$\sigma$ depth')
    plt.legend(scatterpoints=1, loc="upper left", bbox_to_anchor=(1, 1))
    plt.savefig("l45b{0}_cadence.pdf".format(int(b)))

    return times, depths

if __name__ == "__main__":

    outDir = 'LightCurve'
    dbFile = '/Users/ruthangus/Downloads/minion_1016_sqlite.db'
    opsimdb = utils.connectOpsimDb(dbFile)
    resultsDb = db.ResultsDb(outDir=outDir)

    filters = ['u', 'g', 'r', 'i', 'z', 'y']

    # SNR limit (Don't use points below this limit)
    snrLimit = 5.
    # Demand this many points above SNR limit before plotting LC
    nPtsLimit = 6

    # Set RA, Dec for a single point in the sky. in radians.
    ls, bs, ras_deg, decs_deg = np.genfromtxt("coords.txt", skip_header=1).T
    ras_rad = np.radians(ras_deg)
    decs_rad = np.radians(decs_deg)

    for i, ra in enumerate(ras_rad):
        print("b = ", bs[i])
        times, five_sig_depths = get_cadence(ra, decs_rad[i], bs[i], snrLimit,
                                             nPtsLimit, filters, outDir,
                                             opsimdb, resultsDb)
        for j in range(len(filters)):
            data = np.vstack((np.array(times[j]),
                              np.array(five_sig_depths[j])))
            np.savetxt("l45b{0}_{1}_cadence.txt".format(bs[i], filters[j]),
                       data.T)

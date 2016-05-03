import numpy, pylab
import numpy.random
import scipy.interpolate
from SuzPyUtils import filter
from SuzPyUtils.multiplot import *
from SuzPyUtils.norm import *
import glob, atpy, pyfits
import os.path

D2S = 86400
PROT_SUN = 27.0
OMEGA_SUN = 2 * numpy.pi / (27.0 * D2S)

def diffrot_sin2(omega_0, delta_omega, lat):
    return omega_0 - delta_omega * numpy.sin(lat)**2

class spots():
    """Holds parameters for spots on a given star"""
    def __init__(self, dur = None, alpha_med = 0.0001, incl = None, \
                 omega = 2.0, delta_omega = 0.3, diffrot_func = diffrot_sin2, \
                 tau_evol = 5.0, threshold = 0.1,
                 regions = '/Users/aigrain/Soft/idl/diffrot/regions.txt'):
        '''Generate initial parameter set for spots (emergence times
        and initial locations are read from regions file)'''
        # set global stellar parameters which are the same for all spots
        # inclination
        if incl == None:
            self.incl = numpy.arcsin(numpy.random.uniform())
        else:
            self.incl = incl
        # rotation and differential rotation (supplied in solar units)
        self.omega = omega * OMEGA_SUN # in radians
        self.delta_omega = delta_omega * OMEGA_SUN
        self.per_eq = 2 * numpy.pi / self.omega / D2S # in days
        self.per_pole = 2 * numpy.pi / (self.omega - self.delta_omega) / D2S
        self.diffrot_func = diffrot_func
        # spot emergence and decay timescales
        self.tau_em = min(2.0, self.per_eq * tau_evol / 10.0)
        self.tau_decay = self.per_eq * tau_evol
        # read in regions file
        X = numpy.genfromtxt(regions).T
        t0 = X[0].astype(float)
        lat = 0.5*(X[1]+X[3])
        l = lat < 0
	lat = numpy.pi/2. - lat
	lat[l] = -lat[l]
        lon = 0.5*(X[2]+X[4])
        Bem = X[7]
        # keep only spots emerging within specified time-span, with peak B-field > threshold
        if dur == None:
            self.dur = t0.max()
        else:
            self.dur = dur
        l = (t0 < self.dur) * (Bem > threshold)
        self.nspot = l.sum()
        self.t0 = t0[l]
        self.lat = lat[l]
        self.lon = lon[l]
        self.amax = Bem[l] \
            * alpha_med / numpy.median(Bem[l]) # scale to achieve desired median alpha,
                                               # where alpha = spot contrast * spot area

    def calci(self, time, i):
        '''Evolve one spot and calculate its impact on the stellar flux'''
        '''NB: Currently there is no spot drift or shear'''
        # Spot area
        area = numpy.ones(len(time)) * self.amax[i]
        tt = time - self.t0[i]
        l = tt<0
        area[l] *= numpy.exp(-tt[l]**2 / 2. / self.tau_em**2) # emergence
        l = tt>0
        area[l] *= numpy.exp(-tt[l]**2 / 2. / self.tau_decay**2) # decay
        # Rotation rate
        ome = self.diffrot_func(self.omega, self.delta_omega, self.lat[i])
        # Fore-shortening
        phase = ome * time * D2S + self.lon[i]
        beta = numpy.cos(self.incl) * numpy.sin(self.lat[i]) + \
            numpy.sin(self.incl) * numpy.cos(self.lat[i]) * numpy.cos(phase)
        # Differential effect on stellar flux
        dF = - area * beta
        dF[beta < 0] = 0
        return area, ome, beta, dF

    def calc(self, time):
        '''Calculations for all spots'''
        N = len(time)
        M = self.nspot
        area = numpy.zeros((M, N))
        ome = numpy.zeros(M)
        beta = numpy.zeros((M, N))
        dF = numpy.zeros((M, N))
        for i in numpy.arange(M):
            area_i, omega_i, beta_i, dF_i = self.calci(time, i)
            area[i,:] = area_i
            ome[i] = omega_i
            beta[i,:] = beta_i
            dF[i,:] = dF_i
        return area, ome, beta, dF

root = '/Users/aigrain/Data/Kepler/diffrot/'

def mk_lcs():
    dur = 1100
    cad = 5 * 29.4 / 60.0 / 24.0
    X = numpy.genfromtxt('%spar/regions_par.txt' % root).T
    ar = X[1]
    clen = X[2]
    cover = X[3]
    lmin = X[4]
    lmax = X[5]
    rr = X[6]
    nsim = len(ar)
    incl = numpy.arcsin(numpy.sqrt(numpy.random.uniform(0, 1, nsim)))
    n1 = nsim/10
    n2 = nsim-n1
    period = 10.0**(numpy.append(numpy.random.uniform(0, 1, n1), \
                                 numpy.random.uniform(1, 1.7, n2)))
    numpy.random.shuffle(period)
    omega = PROT_SUN / period
    n1 = nsim/3
    n2 = nsim-n1
    delta_omega = numpy.append(10.0**(numpy.random.uniform(-1, 0, n2)), \
                               numpy.zeros(n1))
    numpy.random.shuffle(delta_omega)

    # the following line is wrong - it doesn't actually compute the
    # polar rotation period, but rather a quantity equal to
    # (delta_omega / omega_eq + 1) * p_eq. However, I'm keeping it as
    # it is for the record - the code mk_combined_table.py reads in
    # the spot_par.txt file and corrects this error when making a
    # combined table of all the parameters used to generate every
    # light curve.
    period_pole = period * (1 + delta_omega)

    delta_omega *= omega
    alpha_med = numpy.sqrt(ar) * 3e-4
    tau_evol = 10.0**(numpy.random.uniform(0, 1, nsim))
    f = open('%spar/spot_par.txt' % root, 'w')
    string = '#  N    AR   CLEN  COVER   LMIN   LMAX R   SINI    PEQ   PPOL  A_MED  TAU  NSPOT'
    f.write('%s\n' % string)
    print string
    ee = dofig(1, 1, 3, scale = 1.5)
    for i in numpy.arange(nsim):
        s = spots(incl = incl[i], omega = omega[i], \
                  delta_omega = delta_omega[i], alpha_med = alpha_med[i], \
                  tau_evol = tau_evol[i], regions = '%snoise_free/regions_%04d.txt' % (root, i))
        string = '%04d %5.3f %6.3f %6.3f %6.3f %6.3f %1i %6.2f %6.2f %6.2f %5.2f %5.2f %5i' % \
            (i, ar[i], clen[i], cover[i], lmin[i], lmax[i], rr[i], \
             numpy.sin(incl[i]), period[i], period_pole[i], numpy.log10(alpha_med[i]), \
             tau_evol[i], s.nspot)
        f.write('%s\n' % string)
        print string
        tstart = s.dur - dur
        time = numpy.r_[tstart:s.dur:cad]
        area, ome, beta, dF = s.calc(time)
        X = numpy.genfromtxt('%snoise_free/regions_%04d.txt' % (root, i)).T
        pylab.clf()
        ax1 = doaxes(ee, 1, 3, 0, 0)
        pylab.title('%d AR=%5.3f  CL=%6.3f  R=%1i  sin(i)=%6.2f  PEQ=%6.2f  PPL=%6.2f tau=%5.2f NS=%5i' % \
                    (i, ar[i], clen[i], rr[i], \
                     numpy.sin(incl[i]), period[i], period_pole[i], tau_evol[i], s.nspot))
        for j in numpy.arange(s.nspot):
            pylab.plot(s.t0[j], s.lat[j]*180/numpy.pi, 'ko', markersize = s.amax[j]*(1./3e-4)*5, alpha = 0.5)
        pylab.ylim(-90,90)
        pylab.ylabel('spot lat. (deg)')
        ax2 = doaxes(ee, 1, 3, 0, 1, sharex = ax1)
        pylab.plot(time, area.sum(0), 'k-')
        pylab.ylabel('spot coverage')
        ax3 = doaxes(ee, 1, 3, 0, 2, sharex = ax1)
        pylab.plot(time, dF.sum(0), 'k-')
        pylab.ylabel('delta flux')
        pylab.xlim(time.min(), time.max())
        pylab.xlabel('time (days)')
        pylab.draw()
        pylab.savefig('%snoise_free/lightcurve_%04d.png' % (root, i))
        X = numpy.zeros((2,len(time)))
        X[0,:] = time
        X[1,:] = dF.sum(0)
        numpy.savetxt('%snoise_free/lightcurve_%04d.txt' % (root, i), X.T)
        # raw_input('Next?')
    f.close()
    return

def amps():
    fl = glob.glob('%snoise_free/lightcurve_????.txt' % root)
    n = len(fl)
    amps = numpy.zeros(n)
    for i in numpy.arange(n):
        X = numpy.genfromtxt(fl[i]).T
        dF = X[1]
        amps[i] = dF.max() - dF.min()
        print i, amps[i]
    numpy.savetxt('%spar/amplitudes.txt' % root, amps.T)
    return

def sun():
    dur = 1000
    cad = 29.4 / 60.0 / 24.0
    X = numpy.genfromtxt('%ssun/sun_composite_tsi_20130930.txt' % root).T
    date = X[0]
    t = X[1]
    irr = X[2]
    l = (irr > 0) * (t > 5000)
    t = t[l] - t[l].min()
    date_ref = date[l][0]
    irr = irr[l] / irr[l].max()
    pylab.figure(1)
    pylab.clf()
    pylab.plot(t, irr, 'k-')
    col = ['r','g','b','y','m']
    tstart = [1000, 2100, 2600, 3800, 5000]
    for i in numpy.arange(5):
        time = numpy.r_[tstart[i]:tstart[i]+dur:cad]
        if i == 0:
            y1 = numpy.zeros(len(time)) + 1
            y2 = y1 - 0.0045
        g = scipy.interpolate.interp1d(t, irr, bounds_error = False)
        dF = filter.boxcare(g(time), 10, fill = True)
        pylab.plot(time, dF, c = col[i])
        pylab.fill_between(time, y1, y2, color = col[i], alpha = 0.3)
        X = numpy.zeros((2,len(time)))
        X[0,:] = time - time[0]
        X[1,:] = dF
        numpy.savetxt('%ssun/sun_lightcurve_%02d.txt' % (root, i), X.T)
    pylab.xlim(t.min(), t.max())
    pylab.ylim(0.9955,1)
    pylab.xlabel('Days since %d' % date_ref)
    pylab.ylabel('Flux decrement')
    pylab.savefig('%ssun/sun_lightcurves.png' % root)
    return

def read_kplr(kid, normalise = True, quiet = True):
    files = glob.glob('%s/kplr/kplr%09d-*_llc.fits' % (root, kid))
    nfl = len(files)
    for i in scipy.arange(nfl):
        if quiet == False:
            print files[i]
        hdulist = pyfits.open(files[i])
        # q = hdulist[0].header['QUARTER']
        hdulist.close()
        t = atpy.Table(files[i], type = 'fits', verbose = False)
        tt = t.TIME # + 2454833
        cn = t.CADENCENO
        rf = t.SAP_FLUX
        pf = t.PDCSAP_FLUX
        # qu = scipy.zeros(len(tt), 'int') + q
        if normalise == True:
            l = scipy.isfinite(tt) * scipy.isfinite(rf)
            t0 = tt[l].min()
            p = scipy.polyfit(tt[l] - t0, rf[l], 1)
            rf /= scipy.polyval(p, tt - t0)
            l = scipy.isfinite(tt) * scipy.isfinite(pf)
            t0 = tt[l].min()
            p = scipy.polyfit(tt[l] - t0, pf[l], 1)
            pf /= scipy.polyval(p, tt - t0)
        if i == 0:
            time = tt
            cadence_no = cn
            raw_flux = rf
            pdc_flux = pf
            # quarter = qu
        else:
            time = scipy.append(time, tt)
            cadence_no = scipy.append(cadence_no, cn)
            raw_flux = scipy.append(raw_flux, rf)
            pdc_flux = scipy.append(pdc_flux, pf)
            # quarter = scipy.append(quarter, qu)
    t = atpy.Table()
    t.add_column('time', time, unit = 'BJD')
    t.add_column('cadence_no', cadence_no, dtype = 'int')
    t.add_column('raw_flux', raw_flux, unit = 'e-/s')
    t.add_column('pdc_flux', pdc_flux, unit = 'e-/s')
    # t.add_column('quarter', quarter, dtype = 'int')
    return t

def check_kplr_files():
    kid_lis = numpy.genfromtxt('%skplr/kid.lis' % root)
    files = glob.glob('%skplr/kplr*.fits' % root)
    nfl = len(files)
    kid_downloaded = numpy.zeros(nfl, 'int')
    for i in numpy.arange(nfl):
        ff = os.path.basename(files[i])
        kid_downloaded[i] = int(ff[4:13])
    kid_downloaded = numpy.unique(numpy.sort(kid_downloaded))
    n_missing = len(kid_lis) - len(kid_downloaded)
    l = numpy.in1d(kid_lis, kid_downloaded)
    kid_missing = kid_lis[l == False]
    print len(kid_missing)
    # numpy.savetxt('%skplr/kid_missing_1.lis' % root, kid_missing[:100].T, fmt = '%10d')
    # numpy.savetxt('%skplr/kid_missing_2.lis' % root, kid_missing[100:200].T, fmt = '%10d')
    # numpy.savetxt('%skplr/kid_missing_3.lis' % root, kid_missing[200:300].T, fmt = '%10d')
    # numpy.savetxt('%skplr/kid_missing_4.lis' % root, kid_missing[300:].T, fmt = '%10d')
    return

def inject_noise():
    kid_lis = numpy.genfromtxt('%skplr/kid.lis' % root)
    files = glob.glob('%skplr/kplr*.fits' % root)
    nfl = len(files)
    kid_downloaded = numpy.zeros(nfl, 'int')
    for i in numpy.arange(nfl):
        ff = os.path.basename(files[i])
        kid_downloaded[i] = int(ff[4:13])
    kid_downloaded = numpy.unique(numpy.sort(kid_downloaded))
    n = len(kid_downloaded)
    l = numpy.zeros(n, 'bool')
    for i in numpy.arange(n):
        files = glob.glob('%skplr/kplr%09d*.fits' % (root, kid_downloaded[i]))
        if len(files) >= 16:
            l[i] = True
    kid_downloaded = kid_downloaded[l]
    n = len(kid_downloaded)
    f = open('%spar/matches.txt' % root, 'w')
    string = '#  N KID       TSTART'
    f.write('%s\n' % string)
    print string
    ee = dofig(1, 1, 2, scale = 1.5)
    for i in numpy.arange(n):
        kp = read_kplr(kid_downloaded[i])
        tp = kp.time
        t00 = tp[0]
        tp -= t00
        fp = kp.pdc_flux
        l = scipy.isfinite(tp) * scipy.isfinite(fp)
        tstart = numpy.random.uniform(0,400)
        tend = tstart + 1000
        l *= (tp >= tstart) * (tp < tend)
        tp = tp[l] - tstart
        fp = fp[l]
        X = numpy.genfromtxt('%snoise_free/lightcurve_%04d.txt' % (root, i)).T
        t0 = X[0]
        t01 = t0[0]
        t0 -= t01
        string = '%4d %9d %10.5f' % (i, kid_downloaded[i], tstart+t00-t01)
        f.write('%s\n' % string)
        print string
        f0 = X[1]
        l = scipy.isfinite(t0) * scipy.isfinite(f0)
        t0 = t0[l]
        f0 = f0[l]
        ff = scipy.interpolate.interp1d(t0, f0)(tp)
        fs = ff + fp
        pylab.clf()
        ax1 = doaxes(ee, 1, 2, 0, 0)
        pylab.title('LC %d, KID %d' % (i, kid_downloaded[i]))
        pylab.plot(tp, fp, 'b-')
        pylab.plot(tp, ff+1, 'r-')
        mp,sp = medsig(fp)
        m01 = mymax(ff+1)
        m02 = mymin(ff+1)
        m0r = m01 - m02
        ymax = max(mp + 5 * sp, m01 + 0.1 * m0r)
        ymin = min(mp - 5 * sp, m02 - 0.1 * m0r)
        pylab.ylim(ymin, ymax)
        pylab.ylabel('flux')
        ax3 = doaxes(ee, 1, 2, 0, 1, sharex = ax1)
        pylab.plot(tp, fs, 'k-')
        ms, ss = medsig(fs)
        pylab.ylim(ms - 5 * ss, ms + 5 * ss)
        pylab.ylabel('flux')
        pylab.xlim(tp.min(), tp.max())
        pylab.xlabel('time (days)')
        pylab.draw()
        pylab.savefig('%sfinal/lightcurve_%04d.png' % (root, i))
        X = numpy.zeros((2,len(tp)))
        X[0,:] = tp
        X[1,:] = fs
        numpy.savetxt('%sfinal/lightcurve_%04d.txt' % (root, i), X.T, fmt = '%11.5f')
    f.close()
    return

def transfer_noise_free():
    # kid_lis = numpy.genfromtxt('%skplr/kid.lis' % root)
    # files = glob.glob('%skplr/kplr*.fits' % root)
    # nfl = len(files)
    # kid_downloaded = numpy.zeros(nfl, 'int')
    # for i in numpy.arange(nfl):
    #     ff = os.path.basename(files[i])
    #     kid_downloaded[i] = int(ff[4:13])
    # kid_downloaded = numpy.unique(numpy.sort(kid_downloaded))
    # n = len(kid_downloaded)
    # l = numpy.zeros(n, 'bool')
    # for i in numpy.arange(n):
    #     files = glob.glob('%skplr/kplr%09d*.fits' % (root, kid_downloaded[i]))
    #     if len(files) >= 16:
    #         l[i] = True
    # kid_downloaded = kid_downloaded[l]
    # n = len(kid_downloaded)
    # tp = numpy.r_[0:1000:29.4/60.0/24.0]
    # for j in numpy.arange(1000-n):
    #     i = j+n
    #     print i
    #     X = numpy.genfromtxt('%snoise_free/lightcurve_%04d.txt' % (root, i)).T
    #     t0 = X[0]
    #     t0 -= t0[0]
    #     f0 = X[1]
    #     l = scipy.isfinite(t0) * scipy.isfinite(f0)
    #     t0 = t0[l]
    #     f0 = f0[l]
    #     ff = scipy.interpolate.interp1d(t0, f0)(tp) + 1
    #     X = numpy.zeros((2,len(tp)))
    #     X[0,:] = tp
    #     X[1,:] = ff
    #     numpy.savetxt('%sfinal/lightcurve_%04d.txt' % (root, i), X.T, fmt = '%11.5f')
    for i in numpy.arange(5):
        X = numpy.genfromtxt('%ssun/sun_lightcurve_%02d.txt' % (root, i)).T
        numpy.savetxt('%sfinal/lightcurve_%04d.txt' % (root, i+1000), X.T, fmt = '%11.5f')
    return

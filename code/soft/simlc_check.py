import numpy, pylab

D2S = 86400
PROT_SUN = 27.0
OMEGA_SUN = 2 * numpy.pi / (27.0 * D2S)
root = '/Users/aigrain/Data/Kepler/diffrot/'

X = numpy.genfromtxt('%snoise_free/regions_par.txt' % root).T
lmin = X[4]
lmax = X[5]
nsim = len(lmin)
incl = numpy.arcsin(numpy.sqrt(numpy.random.uniform(0, 1, nsim)))
n1 = nsim/10
n2 = nsim-n1
period_eq = 10.0**(numpy.append(numpy.random.uniform(0, 1, n1), \
                                numpy.random.uniform(1, 1.7, n2)))
numpy.random.shuffle(period_eq) # in days
omega_eq_solar = PROT_SUN / period_eq # relative to solar rotation rate
n1 = nsim/3
n2 = nsim-n1
delta_per_rel = numpy.append(10.0**(numpy.random.uniform(-1, 0, n2)), \
                             numpy.zeros(n1))
numpy.random.shuffle(delta_per_rel) # relative
period_pole = period_eq * (1 + delta_per_rel)
delta_omega_solar = delta_per_rel * omega_eq_solar # relative to solar rotation rate
omega_eq = omega_eq_solar * OMEGA_SUN # in radians
delta_omega = delta_omega_solar * OMEGA_SUN # in radians
omega_pole = omega_eq - delta_omega # in radians
period_eq_2 = 2 * numpy.pi / omega_eq / D2S # in days
period_pole_2 = 2 * numpy.pi / omega_pole / D2S # in days

x = period_pole_2
y = period_eq / (1 - period_eq * delta_per_rel / PROT_SUN)
pylab.clf()
pylab.plot(x, y, 'k.')
pylab.plot([0,200],[0,200])
pylab.xlim(0,200)
pylab.ylim(0,200)

for i in numpy.arange(nsim):
    print 'Saved in output file:'
    print 'Period_eq: ', period_eq[i]
    print 'Period_pole: ', period_pole[i]
    print 'Differential rotation in terms of relative period'
    print 'Delta_per_rel: ', delta_per_rel[i]
    print 'Rotation rate in solar units'
    print 'Omega_eq_solar: ', omega_eq_solar[i]
    print 'Delta_omega_solar: ', delta_omega_solar[i]
    print 'Rotation rate in radians'
    print 'Omega_eq: ', omega_eq[i]
    print 'Delta_omega: ', delta_omega[i]
    print 'Omega_pole: ', omega_pole[i]
    print 'Computed from Omega_eq and Omega_pole'
    print 'Period_eq: ', period_eq_2[i]
    print 'Period_pole: ', period_pole_2[i]
    print ''
    raw_input('Next?')


import numpy
import pylab
import atpy

root = '/Users/aigrain/Data/Kepler/diffrot'

# read table of spot parameters
par = atpy.Table('%s/par/spot_par.txt' % root, type = 'ascii')

# compute delta_omega_rel and related parameters
peq = par.PEQ # this is the equatorial period in days
ppol_sav = par.PPOL # this is not the actual polar period but can be
                    # used to find delta_per_rel
delta_omega_rel = ppol_sav / peq - 1.0  # this is in fact the relative
                                        # difference between the
                                        # equatorial and polar rotation
                                        # ***rates***, delta_omega / omega_eq
omega_eq = 2 * numpy.pi / peq # equatorial rotation rate in radians / day
delta_omega = delta_omega_rel * omega_eq # absolute differential rotation, in radians / day
omega_pole = omega_eq - delta_omega # polar rotation rate in radians / day
ppol = 2 * numpy.pi / omega_pole # polar rotation period in days
l_min = par.LMIN * numpy.pi / 180.0 # min latitude of spots, in radians
l_max = par.LMAX * numpy.pi / 180.0 # max latitude of spots, in radians
omega_max =  omega_eq - delta_omega * numpy.sin(l_min)**2 # max/min observable rotation rate, ...
omega_min =  omega_eq - delta_omega * numpy.sin(l_max)**2 # ... in radians / day
p_max = 2 * numpy.pi / omega_min # min / max observable rotation period, ...
p_min = 2 * numpy.pi / omega_max # ... in days

# update table
par.add_column('DELTA_OMEGA', delta_omega)
par.add_column('DELTA_OMEGA_REL', delta_omega_rel)
par.add_column('OMEGA_EQ', omega_eq)
par.add_column('OMEGA_POL', omega_pole)
par.PPOL = ppol
par.add_column('OMEGA_MIN', omega_min)
par.add_column('OMEGA_MAX', omega_max)
par.add_column('P_MIN', p_min)
par.add_column('P_MAX', p_max)

# read amplitudes
amps = numpy.genfromtxt('%s/par/amplitudes.txt' % root).flatten()
par.add_column('AMP', amps)

# read amplitudes
matches = atpy.Table('%s/par/matches.txt' % root, type = 'ascii')
kid = matches.KID[:len(par)]
kid[750:] = 0
par.add_column('KID', kid)

par.write('%s/par/final_table.txt' % root, type = 'ascii', overwrite = True)

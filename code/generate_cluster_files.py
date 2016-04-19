from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from teff_bv import teff2bv

# NGC 5316: 90 G dwarfs, 240 K dwarfs, 1650 M dwarfs
# Age = .170 Gyrs, distance = 1208 pc, dm = 10.410
# NGC 2477: 90 G dwarfs, 240 K dwarfs, 1650 M dwarfs
# Age = .822 Gyrs, distance = 1450, dm = 10.807
# IC 4651: 75 G dwarfs, 1500 M dwarfs.
# Age = 1.778 Gyrs, distance = 888, dm = 9.742

# G = 5200, 6000
# K = 3700, 5200
# M = 2400, 3700

# data = np.vstack((logAges, bvs, logTeff, rmag))

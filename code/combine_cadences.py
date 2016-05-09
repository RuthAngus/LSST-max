import numpy as np
import matplotlib.pyplot as plt

def collate(b):

    filters = ["u", "g", "r", "i", "z", "y"]
    all_times, all_depths = [], []
    for i, f in enumerate(filters):
        times, depths = \
            np.genfromtxt("l45b{0}_{1}_cadence.txt".format(b, f)).T
        all_times.append(times)
        all_depths.append(depths)
    t = np.array([i for j in all_times for i in j])
    d = np.array([i for j in all_depths for i in j])
    inds = np.argsort(t)
    t, d = t[inds], d[inds]
    data = np.vstack((t, d))
    np.savetxt("b{0}_cadence.txt".format(b), data.T)

if __name__ == "__main__":

    bs = [-10, -20, -40, -60, -80]
    for b in bs:
        collate(float(b))

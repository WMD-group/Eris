# Script for scatter contour plot based on: https://stackoverflow.com/questions/18764814/make-contour-of-scatter

import matplotlib.pyplot as plt
import numpy as np
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata

def grid(x, y, z, resX=100, resY=100):
    "Convert 3 column data to matplotlib grid"
    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z = plt.mlab.griddata(x, y, z, xi, yi, interp='linear')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

# Output file from Eris to plot (can later make loop over T and savefig for all plots)
data_file = "PotentialCuZn_T_0000K_slice_z=07.dat"

# Reading in data from eris output file
CuZnSlice = np.genfromtxt(data_file, delimiter = ' ')
x_vals = CuZnSlice[:,0]
y_vals = CuZnSlice[:,1]
pots = CuZnSlice[:,2] 

X, Y, Z = grid(x_vals, y_vals, pots)
plt.contourf(X, Y, Z)

plt.xlabel('X (lattice units)')
plt.ylabel('Y (lattice units)')
plt.colorbar()

plt.show()
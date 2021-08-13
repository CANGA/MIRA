#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 11:42:40 2019

Tests the computation of gradients on target meshes. Assumes that data is loaded
into the workspace by running CANGAMetricsDriver.py CHECKS ON RLL1deg.

@author: jeguerra
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Compute gradients here if necessary...
nshape = (180, 360)

# Compute the magnitude of the gradients
gradMag1 = np.linalg.norm(gradientsOnTM[0], axis=0)
gradMag2 = np.linalg.norm(gradientsOnTM[1], axis=0)

# Fetch the gradients
varpST = np.reshape(varST, nshape)
gradST = np.reshape(gradMag1, nshape)
gradS2T = np.reshape(gradMag2, nshape)

# Make 1D arrays for tri-surface plots
x = np.reshape(cellCoordT[0,:], nshape)
y = np.reshape(cellCoordT[1,:], nshape)
#z = cellCoordT[2,:]

#X, Y = np.meshgrid(x, y)

# Make a nice 3D plot of the variable
fig = plt.figure()
cs = plt.contourf(varpST, 200, cmap='RdGy')
plt.colorbar()
#plt.colorbar(m, boundaries=np.arange(0,3.1,.5))

# Make a nice 3D plot of the gradients
vmin = 0.0
vmax = 5000.0
fig = plt.figure()
#cs = plt.contourf(gradS2T, 200, cmap='RdGy', vmin = vmin, vmax = vmax)
cs = plt.contourf(gradS2T, 100, cmap='RdGy')
m = plt.cm.ScalarMappable(cmap='RdGy')
m.set_array(gradS2T)
m.set_clim(vmin, vmax)
plt.colorbar(m)
#plt.colorbar(m, boundaries=np.arange(0,3.1,.5))
plt.show()




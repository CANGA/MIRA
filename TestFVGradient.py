#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 11:42:40 2019

Tests the computation of gradients on target meshes. Assumes that data is loaded
into the workspace by running CANGAMetricsDriver.py CHECKS ON RLL.

@author: jeguerra
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Compute gradients here if necessary...
nshape = (180, 360)

# Fetch the gradients
varpST = np.reshape(varST, nshape)
gradST = np.reshape(gradientsOnTM[0], nshape)
gradS2T = np.reshape(gradientsOnTM[1], nshape)

# Make 1D arrays for tri-surface plots
x = np.reshape(cellCoordT[0,:], nshape)
y = np.reshape(cellCoordT[1,:], nshape)
#z = cellCoordT[2,:]

#X, Y = np.meshgrid(x, y)

# Make a nice 3D plot of the variable
fig = plt.figure()
cs = plt.contourf(varpST, 100, cmap='RdGy')
plt.colorbar()
#plt.colorbar(m, boundaries=np.arange(0,3.1,.5))

# Make a nice 3D plot of the gradients
vmin = 50.0
vmax = 500.0
fig = plt.figure()
cs = plt.contourf(gradS2T, 100, cmap='RdGy', vmin = vmin, vmax = vmax)
m = plt.cm.ScalarMappable(cmap='RdGy')
m.set_array(gradS2T)
m.set_clim(vmin, vmax)
plt.colorbar(m)
#plt.colorbar(m, boundaries=np.arange(0,3.1,.5))
plt.show()




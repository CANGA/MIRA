#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 10:33:59 2019

@author: TempestGuerra
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ASSUMES THE COORDINATES ARE ALREADY IN THE WORKSPACE
ax.scatter(varCoordT[0,:], varCoordT[1,:], varCoordT[2,:])
plt.show()
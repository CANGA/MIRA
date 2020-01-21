#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 10:48:07 2020

Checks modal convergence of SH reconstruction for Terrain field

@author: TempestGuerra
"""
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

# Bring in the results
f32id = Dataset('/Users/TempestGuerra/Desktop/Remapping-Intercomparison/testdata_NM32outCSne60_TPO.nc')
f64id = Dataset('/Users/TempestGuerra/Desktop/Remapping-Intercomparison/testdata_NM64outCSne60_TPO.nc')
f128id = Dataset('/Users/TempestGuerra/Desktop/Remapping-Intercomparison/testdata_NM128outCSne60_TPO.nc')
f256id = Dataset('/Users/TempestGuerra/Desktop/Remapping-Intercomparison/testdata_NM256outCSne60_TPO.nc')
f512id = Dataset('/Users/TempestGuerra/Desktop/Remapping-Intercomparison/testdata_NM512outCSne60_TPO.nc')

# Compute absolute differences
ref512 = f512id.variables['Topography'][:]
dtp_32 = np.abs(f32id.variables['Topography'][:] - ref512)
dtp_64 = np.abs(f64id.variables['Topography'][:] - ref512)
dtp_128 = np.abs(f128id.variables['Topography'][:] - ref512)
dtp_256 = np.abs(f256id.variables['Topography'][:] - ref512)

# Get norms and plot errors
errors = [np.linalg.norm(dtp_32) / len(dtp_32), \
          np.linalg.norm(dtp_64) / len(dtp_64), \
          np.linalg.norm(dtp_128) / len(dtp_128), \
          np.linalg.norm(dtp_256) / len(dtp_256)]
nmodes = [32, 64, 128, 256]

plt.plot(nmodes, errors); plt.yscale('log')
#plt.title('Terrain Field by SH Convergence: Centroid Sampling')
#plt.title('Terrain Field by SH Convergence: 2nd Order Sampling')
plt.title('Terrain Field by SH Convergence: 4th Order Sampling')
plt.xlabel('Number of SH Modes')
plt.ylabel('L2 Error WRT N=512')
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:45:49 2019

RUN THE FOLLOWING FIRST:
       
python3 CANGAMetricsDriver.py -v Psi 
--ss testdata_CSne16_np4_1.nc 
--st testdata_ICO16_np4_1.nc 
--s2t testdata_ICO16_np4_1.nc 
--sm meshes/outCSne16.g_RIPPED 
--tm meshes/outICO16.g_RIPPED 
--smc 1 --tmc 1

Data in "testdata" files is not used. This run is to load the mesh data.

@author: jeguerra
"""
import pyshtools
import numpy as np
import math as mt
from computeAreaIntegral import computeAreaIntegral

'''
#%% EVALUATE TEST = 1 DATA AND EXPAND IN SH
NX = 720
NY = 360
latgrid, longrid = np.meshgrid(np.linspace(-mt.pi, mt.pi, NY), \
                               np.linspace(0.0, 2.0 * mt.pi, NX), indexing='ij')

test1 = np.zeros((NY,NX))
test2 = np.zeros((NY,NX))
for ii in range(NY):
       for jj in range(NX):
              test1[ii,jj] = (2.0 + mt.cos(latgrid[ii,jj]) * \
                            mt.cos(latgrid[ii,jj]) * \
                            mt.cos(2.0 * longrid[ii,jj])) # TEST 1
                     
              test2[ii,jj] = (2.0 + (np.sin(2.0 * latgrid[ii,jj]))**16.0 * \
                              np.cos(16.0 * longrid[ii,jj])) # TEST 2
              
# Expand CS16 in SH
coeffs = pyshtools.expand.SHExpandDH(test2, sampling=2)
clm = pyshtools.SHCoeffs.from_array(coeffs)
'''
AREAFV_S = 0.0
AREAFV_T = 0.0
order = 10

NC = len(varConS)
areaFV = np.zeros((NC,1))
for ii in range(NC):
       # Compute areas by triangular quadrature on source FV mesh
       cdex = varConS[ii,:] - 1
       thisCell = varCoordS[:,cdex.astype(int)]
       areaFV[ii] = computeAreaIntegral(None, thisCell, order, False, False)
       AREAFV_S += areaFV[ii]
       
NC = len(varConT)
areaFV = np.zeros((NC,1))
for ii in range(NC):
       # Compute areas by triangular quadrature on target FV mesh
       cdex = varConT[ii,:] - 1
       thisCell = varCoordT[:,cdex.astype(int)]
       areaFV[ii] = computeAreaIntegral(None, thisCell, order, False, False)
       AREAFV_T += areaFV[ii]

#%%
#print("Unit Sphere Area, Built-in reference value: %16.15e" %(4.0 * mt.pi))
print("Test 2 Sphere Integral, computed on CS16 @ " + str(order) + "th order: %16.15e" %AREAFV_S)
print("Test 2 Sphere Integral, computed on ICO16 @ " + str(order) + "th order: %16.15e" %AREAFV_T)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TEST FOR QUADRATURE CONVERGENCE
@author: jeguerra
"""
import pyshtools
import numpy as np
import math as mt
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from computeAreaIntegral import computeAreaIntegral
import computeSphericalCartesianTransforms as sphcrt

#%% USER INPUT! ** CHANGE THESE TO RUN **

# Flag to use SH expansions (True) or evaluate directly (False)
USESH = True
meshesFolder = '/Users/TempestGuerra/Desktop/Remapping-Intercomparison/convergenceMeshes/'

#%% ANALYTICAL TEST FUNCTIONALS (PASSED INTO AREA INTEGRAL FUNCTION)

def test1(lonRAD, latRAD):
       val = (2.0 + mt.cos(latRAD * \
              mt.cos(latRAD) * \
              mt.cos(2.0 * lonRAD)))
              
       return val

def test2(lonRAD, latRAD):
       val = (2.0 + (np.sin(2.0 * latRAD))**16.0 * \
               np.cos(16.0 * lonRAD))
              
       return val

#%% EVALUATE THE TWO SMOOTH FIELDS AND GENERATE SH EXPANSIONS
if USESH:
       NX = 720
       NY = 360
       latgrid, longrid = np.meshgrid(np.linspace(-mt.pi, mt.pi, NY), \
                                      np.linspace(0.0, 2.0 * mt.pi, NX), indexing='ij')
       
       valTest1 = np.zeros((NY,NX))
       valTest2 = np.zeros((NY,NX))
       for ii in range(NY):
              for jj in range(NX):
                     valTest1[ii,jj] = (2.0 + mt.cos(latgrid[ii,jj]) * \
                                   mt.cos(latgrid[ii,jj]) * \
                                   mt.cos(2.0 * longrid[ii,jj])) # TEST 1
                            
                     valTest2[ii,jj] = (2.0 + (np.sin(2.0 * latgrid[ii,jj]))**16.0 * \
                                     np.cos(16.0 * longrid[ii,jj])) # TEST 2
                     
       # Get the SH expansion object for each of the tests (clm1 and clm2)
       coeffs = pyshtools.expand.SHExpandDH(valTest1, sampling=2)
       clm1 = pyshtools.SHCoeffs.from_array(coeffs)
       coeffs = pyshtools.expand.SHExpandDH(valTest2, sampling=2)
       clm2 = pyshtools.SHCoeffs.from_array(coeffs)

#%% READ IN THE MESH FILES AND GET COORDINATE DATA

meshRes = ('16','32','64','128')
meshCoords = []
meshCells = []
for ff in range(4):
       thisMeshFile = meshesFolder + 'outCSne' + meshRes[ff] + '.g'
       thisFid = Dataset(thisMeshFile,'r')
       meshCoords.append(thisFid.variables['coord'][:])
       meshCells.append(thisFid.variables['connect1'][:])
       thisFid.close()
       print('Mesh information for CSne' + meshRes[ff] + '...DONE!')

#%% COMPUTE CENTROID VALUE AND CELL AVERAGES (6th order quadrature) FOR EACH MESH

differences = []
NC = []
order = 4

for mm in range(4):
       # Get mesh information at this resolution
       thisCoord = meshCoords[mm]
       thisCell = meshCells[mm]
       thisCentroids = sphcrt.computeCentroids(thisCell, thisCoord)
       thisCentroidsLL = sphcrt.computeCart2LL(thisCentroids)
       thisCentroidsLL_RAD = mt.pi / 180.0 * thisCentroidsLL
       numCells = len(thisCell)
       NC.append(numCells)
       
       # Get the centroid values directly from SH expansion
       if USESH:
              thisValues1 = clm1.expand(lon=thisCentroidsLL[:,0], lat=thisCentroidsLL[:,1])
              thisValues2 = clm2.expand(lon=thisCentroidsLL[:,0], lat=thisCentroidsLL[:,1])
       else:
              thisValues1 = np.zeros(numCells)
              thisValues2 = np.zeros(numCells)
              for cc in range(numCells):
                     thisValues1[cc] = test1(thisCentroidsLL_RAD[cc,0], thisCentroidsLL_RAD[cc,1])
                     thisValues2[cc] = test2(thisCentroidsLL_RAD[cc,0], thisCentroidsLL_RAD[cc,1])
       
       # Get the area cell averages for the functions
       thisAreas1 = np.zeros(numCells)
       thisAreas2 = np.zeros(numCells)
       for cc in range(numCells):
              cdex = thisCell[cc,:] - 1
              aCell = thisCoord[:,cdex.astype(int)]
              
              # Get the area average of each function
              if USESH:
                     thisAreas1[cc] = computeAreaIntegral(clm1, aCell, order, True, False)
                     thisAreas2[cc] = computeAreaIntegral(clm2, aCell, order, True, False)
              else:
                     thisAreas1[cc] = computeAreaIntegral(test1, aCell, order, True, False)
                     thisAreas2[cc] = computeAreaIntegral(test2, aCell, order, True, False)
              
       # Now store the differences in tuples
       differences.append((np.abs(thisAreas1 - thisValues1), \
                           np.abs(thisAreas2 - thisValues2)))
       print('Computed differences for CSne' + meshRes[mm] + '...DONE!')
       
#%% COMPUTE NORMS OF DIFFERENCES AND PLOT

ii = 0
Linf_test1 = []
Linf_test2 = []
for res in differences:
       # Compute L_inf norm for each test on each mesh
       Linf_test1.append((1.0 / NC[ii]) * np.linalg.norm(res[0], ord=np.inf))
       Linf_test2.append((1.0 / NC[ii]) * np.linalg.norm(res[1], ord=np.inf))
       
       print(meshRes[ii], Linf_test1[ii], Linf_test2[ii])
       ii += 1

#%% MAKE THE CONVERGENCE PLOTS (FINGERS CROSSED...)
meshDelta = np.array([1.0, 0.5, 0.25, 0.125])
scale = 0.5 * (Linf_test1[0] + Linf_test2[0])
order2 = scale * np.power(meshDelta,2.0)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, Linf_test1, 'bs-')
plt.plot(meshDelta, Linf_test2, 'r+-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Test 1', 'Test 2'), loc='upper right')
plt.title('Grid Convergence Test: Centroid Sample to 4th Order Quadrature')
plt.xlabel('Normalized Grid Spacing')
plt.ylabel('Infinity Norm Error')
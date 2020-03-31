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
from computeStandardNorms import computeStandardNorms

#%% USER INPUT! ** CHANGE THESE TO RUN **

# Flag to set norm computation
byStandardNorms = False

# Flag to use SH expansions (True) or evaluate directly (False)
USESH = False
# RUN: meshes/MakeConvergenceMeshesTempestRemap.sh (Change location of tempestremap executables)
# Change the following variable to the path of the meshes created above
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
values = []
references = []
areas = []
NC = []
order = 4

for mm in range(4):
       # Get mesh information at this resolution
       thisCoord = meshCoords[mm]
       thisCell = meshCells[mm]
       thisCentroids = sphcrt.computeCentroids(thisCell, thisCoord)
       thisCentroidsLL_RAD = sphcrt.computeCart2LL(thisCentroids)
       thisCentroidsLL_DEG = 180.0 / mt.pi * thisCentroidsLL_RAD
       numCells = len(thisCell)
       NC.append(numCells)
       
       # Get the centroid values directly from SH expansion or NOT
       if USESH:
              thisValues1 = clm1.expand(lon=thisCentroidsLL_DEG[:,0], lat=thisCentroidsLL_DEG[:,1])
              thisValues2 = clm2.expand(lon=thisCentroidsLL_DEG[:,0], lat=thisCentroidsLL_DEG[:,1])
       else:
              thisValues1 = np.zeros(numCells)
              thisValues2 = np.zeros(numCells)
              for cc in range(numCells):
                     thisValues1[cc] = test1(thisCentroidsLL_RAD[cc,0], thisCentroidsLL_RAD[cc,1])
                     thisValues2[cc] = test2(thisCentroidsLL_RAD[cc,0], thisCentroidsLL_RAD[cc,1])
       
       # Get the area cell averages for the functions
       thisAvgs1 = np.zeros(numCells)
       thisAvgs2 = np.zeros(numCells)
       thisAreas = np.zeros(numCells)
       for cc in range(numCells):
              cdex = thisCell[cc,:] - 1
              aCell = thisCoord[:,cdex.astype(int)]
              
              # Get the area average of each function
              if USESH:
                     thisAvgs1[cc] = computeAreaIntegral(clm1, aCell, order, True, False)
                     thisAvgs2[cc] = computeAreaIntegral(clm2, aCell, order, True, False)
              else:
                     thisAvgs1[cc] = computeAreaIntegral(test1, aCell, order, True, False)
                     thisAvgs2[cc] = computeAreaIntegral(test2, aCell, order, True, False)
              
              # Compute the areas
              thisAreas[cc] = computeAreaIntegral(None, aCell, order, False, True)
       
       # Store the differences in tuples
       differences.append((np.abs(thisAvgs1 - thisValues1), \
                           np.abs(thisAvgs2 - thisValues2)))
              
       #Store the low order evaluation
       values.append((thisValues1, thisValues2))
              
       # Store the reference (high order cell averages)
       references.append((thisAvgs1, thisAvgs2))
       
       # Store the areas
       areas.append(thisAreas)
       
       print('Computed differences for CSne' + meshRes[mm] + '...DONE!')
       
#%% COMPUTE STANDARD NORMS
if byStandardNorms:
       ii = 0
       L1_test1 = []
       L1_test2 = []
       L2_test1 = []
       L2_test2 = []
       Li_test1 = []
       Li_test2 = []
       
       for diff in differences:
              # Compute standard norms for tests 1 and 2
              L1_t1, L2_t1, Li_t1 = computeStandardNorms(meshCells[ii], values[ii][0], references[ii][0], areas[ii], None, False)
              L1_t2, L2_t2, Li_t2 = computeStandardNorms(meshCells[ii], values[ii][1], references[ii][1], areas[ii], None, False)
              
              L1_test1.append(L1_t1)
              L1_test2.append(L1_t2)
              
              L2_test1.append(L2_t1)
              L2_test2.append(L2_t2)
              
              Li_test1.append(Li_t1)
              Li_test2.append(Li_t2)
              
              print('L1 Norms by computeStandardNorms()...')
              print(meshRes[ii], L1_test1[ii], L1_test2[ii])
              ii += 1
else:
       #% COMPUTE STANDARD NORMS OF DIFFERENCES (BY INDEPENDENT MEANS...)
       ii = 0
       L1_test1 = []
       L1_test2 = []
       L2_test1 = []
       L2_test2 = []
       Li_test1 = []
       Li_test2 = []
       
       for diff in differences:
              
              # Compute L1 Norm
              L1_test1.append(diff[0].dot(areas[ii]) / (references[ii][0]).dot(areas[ii]))
              L1_test2.append(diff[1].dot(areas[ii]) / (references[ii][1]).dot(areas[ii]))
              
              # Compute L2 Norm
              L2_test1.append(np.power(diff[0],2.0).dot(areas[ii]) / (np.power(references[ii][0],2.0)).dot(areas[ii]))
              L2_test2.append(np.power(diff[1],2.0).dot(areas[ii]) / (np.power(references[ii][1],2.0)).dot(areas[ii]))
              
              # Compute max Norm
              Li_test1.append(np.amax(diff[0] * areas[ii]) / np.amax(references[ii][0] * areas[ii]))
              Li_test2.append(np.amax(diff[1] * areas[ii]) / np.amax(references[ii][1] * areas[ii]))
              
              print('L1 Norms by independent code...')
              print(meshRes[ii], L1_test1[ii], L1_test2[ii])
              ii += 1

#%% MAKE THE CONVERGENCE PLOTS (FINGERS CROSSED...)
meshDelta = np.array([1.0, 0.5, 0.25, 0.125])
plt.figure(figsize=(6.0, 10.0))

# L1 Norm plot
scale = 0.5 * (L1_test1[0] + L1_test2[0])
order2 = scale * np.power(meshDelta,2.0)
plt.subplot(3,1,1)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, L1_test1, 'bs-')
plt.plot(meshDelta, L1_test2, 'r+-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Test 1', 'Test 2'), loc='upper right')
plt.title('Grid Convergence Test: Centroid Sample to 4th Order Quadrature')
plt.ylabel('L1 Norm Error')

# L2 Norm plot
scale = 0.5 * (L2_test1[0] + L2_test2[0])
order2 = scale * np.power(meshDelta,2.0)
plt.subplot(3,1,2)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, L2_test1, 'bs-')
plt.plot(meshDelta, L2_test2, 'r+-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Test 1', 'Test 2'), loc='upper right')
plt.ylabel('L2 Norm Error')

# Li Norm plot
scale = 0.5 * (Li_test1[0] + Li_test2[0])
order2 = scale * np.power(meshDelta,2.0)
plt.subplot(3,1,3)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, Li_test1, 'bs-')
plt.plot(meshDelta, Li_test2, 'r+-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Test 1', 'Test 2'), loc='upper right')
plt.xlabel('Normalized Grid Spacing')
plt.ylabel('Linf Norm Error')
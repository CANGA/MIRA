#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TEST FOR QUADRATURE CONVERGENCE FROM SAMPLED Topography DATA
@author: jeguerra
"""
import numpy as np
import math as mt
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from computeAreaIntegral import computeAreaIntegral
import computeSphericalCartesianTransforms as sphcrt
from computeStandardNorms import computeStandardNorms

#%% USER INPUT! ** CHANGE THESE TO RUN **
isEnhanced = True

# RUN: meshes/MakeConvergenceMeshesTempestRemap.sh (Change location of tempestremap executables)
# Change the following variable to the path of the meshes created above
meshesFolder = '/Users/TempestGuerra/Desktop/Remapping-Intercomparison/convergenceMeshes/'

# RUN: ./processConvergenceMeshes.sh
# RUN: ./makeGridConvergenceData_N64_TPO.sh

#%% READ IN THE MESH FILES AND GET COORDINATE DATA

meshRes = ('16','32','64','128')
meshCoords = []
meshCells = []
areas = []
stencils = []
for ff in range(4):
       thisMeshFile = meshesFolder + 'outCSne' + meshRes[ff] + '_enhanced.g'
       thisFid = Dataset(thisMeshFile,'r')
       meshCoords.append(thisFid.variables['coord'][:])
       meshCells.append(thisFid.variables['connect1'][:])
       
       thisCoord = meshCoords[ff]
       thisCell = meshCells[ff]
       numCells = len(thisCell)
       
       if isEnhanced:
              areas.append(thisFid.variables['cell_area'][:])
              stencils.append(thisFid.variables['cell_edge_adjacency'][:])
       else:
              print('MUST RUN WITH PREPROCESSED MESH FILES FOR ADJACENCY DATA!')
              stencils.append(None)
              
              thisAreas = np.zeros(numCells)
              for cc in range(numCells):
                     cdex = thisCell[cc,:] - 1
                     aCell = thisCoord[:,cdex.astype(int)]
                     # Compute the areas
                     thisAreas[cc] = computeAreaIntegral(None, aCell, 6, False, True)
                     
              # Store the areas
              areas.append(thisAreas)
       
       thisFid.close()
       print('Mesh information for CSne' + meshRes[ff] + '...DONE!')

#%% COMPUTE CENTROID VALUE AND CELL AVERAGES (6th order quadrature) FOR EACH MESH

differences = []
values = []
references = []
NC = []
order = 4

for mm in range(4):
       # Get mesh information at this resolution
       thisCoord = meshCoords[mm]
       thisCell = meshCells[mm]
       thisCentroids = sphcrt.computeCentroids(thisCell, thisCoord)
       thisCentroidsLL_RAD = sphcrt.computeCart2LL(thisCentroids)
       thisCentroidsLL_DEG = 180.0 / mt.pi * thisCentroidsLL_RAD
       
       # Open the sample files for reading
       fnameO1 = 'sample_NM64_O1_outCSne' + str(meshRes[mm]) + '_A2.nc'
       fidO1 = Dataset(fnameO1, 'r')
       fnameO4 = 'sample_NM64_O8_outCSne' + str(meshRes[mm]) + '_A2.nc'
       fidO4 = Dataset(fnameO4, 'r')
       
       # Get the sampled data
       thisValues = fidO1.variables['AnalyticalFun2'][:]
       thisAvgs = fidO4.variables['AnalyticalFun2'][:]
       
       # Store the differences
       differences.append(np.abs(thisAvgs - thisValues))
              
       #Store the low order evaluation
       values.append(thisValues)
              
       # Store the reference (high order cell averages)
       references.append(thisAvgs)
              
       print('Computed differences for CSne' + meshRes[mm] + '...DONE!')
       
#%% COMPUTE LOCAL GRADIENTS ON BOTH TYPES OF SAMPLED DATA
grad_values = []
grad_references = []
grad_differences = []
from computeGradientFV3 import computeGradientFV3
for mm in range(4):
       # Get gradients of centroid sampled data
       thisGradValues = computeGradientFV3(values[mm], meshCells[mm], meshCoords[mm], stencils[mm])
       # Get gradients of cell-averaged sampled data
       thisGradAvgs = computeGradientFV3(references[mm], meshCells[mm], meshCoords[mm], stencils[mm])
       
       # Store the differences
       grad_differences.append(np.abs(thisGradAvgs - thisGradValues))
              
       #Store the low order evaluation
       grad_values.append(thisGradValues)
              
       # Store the reference (high order cell averages)
       grad_references.append(thisGradAvgs)
              
       print('Computed gradient differences for CSne' + meshRes[mm] + '...DONE!')
       
#%% COMPUTE STANDARD NORMS
ii = 0
L1_test1 = []
L2_test1 = []
Li_test1 = []

for diff in differences:
       # Compute standard norms for tests 1 and 2
       L1_t1, L2_t1, Li_t1 = computeStandardNorms(meshCells[ii], values[ii], references[ii], areas[ii], None, False)
       
       L1_test1.append(L1_t1)
       L2_test1.append(L2_t1)
       Li_test1.append(Li_t1)
       
       print('L1 Norms by computeStandardNorms()...')
       print(meshRes[ii], L1_test1[ii])
       ii += 1
       
#%% COMPUTE GRADIENT PRESERVATION METRICS
ii = 0
H1_test1 = []
H1_2_test1 = []

from computeGradientPreserveMetrics import computeGradientPreserveMetrics
for diff in grad_differences:
       varsOnTM = [references[ii], values[ii]]
       gradientsOnTM = [grad_references[ii], grad_values[ii]]

       # Gradient preservation
       H1, H1_2 = computeGradientPreserveMetrics(meshCells[ii], gradientsOnTM, varsOnTM, areas[ii], areas[ii], False)
       
       H1_test1.append(H1)
       H1_2_test1.append(H1_2)
       ii += 1

#%% MAKE THE CONVERGENCE PLOTS (STANDARD NORMS)
meshDelta = np.array([1.0, 0.5, 0.25, 0.125])
plt.figure(figsize=(6.0, 10.0))

# L1 Norm plot
scale = L1_test1[0]
order2 = scale * np.power(meshDelta,2.0)
plt.subplot(3,1,1)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, L1_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Functional',), loc='upper right')
plt.title('Grid Convergence: Centroid Sample to 4th Order Quadrature')
plt.ylabel('L1 Norm Error')

# L2 Norm plot
scale = L2_test1[0]
order2 = scale * np.power(meshDelta,2.0)
plt.subplot(3,1,2)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, L2_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Functional'), loc='upper right')
plt.ylabel('L2 Norm Error')

# Li Norm plot
scale = Li_test1[0]
order2 = scale * np.power(meshDelta,2.0)
plt.subplot(3,1,3)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, Li_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('2nd Order', 'Functional'), loc='upper right')
plt.xlabel('Normalized Grid Spacing')
plt.ylabel('Linf Norm Error')

#%% MAKE THE CONVERGENCE PLOTS (GRADIENT NORMS)
meshDelta = np.array([1.0, 0.5, 0.25, 0.125])
plt.figure(figsize=(6.0, 10.0))

# H1 Norm plot
scale = H1_test1[0]
order2 = scale * np.power(meshDelta,1.5)
plt.subplot(2,1,1)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, H1_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('1.5th Order', 'Functional',), loc='upper right')
plt.title('Grid Convergence: Centroid Sample to 4th Order Quadrature')
plt.ylabel('H1 Gradient Norm Error')

# H1_2 Norm plot
scale = H1_2_test1[0]
order2 = scale * np.power(meshDelta,1.5)
plt.subplot(2,1,2)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, H1_2_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(b=None, which='major', axis='both', color='k', linestyle='--', linewidth=0.5)
plt.legend(('1.5th Order', 'Gradient'), loc='upper right')
plt.ylabel('H-1/2 Gradient Semi-Norm Error')
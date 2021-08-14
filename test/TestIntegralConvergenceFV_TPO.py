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
import computeSphericalCartesianTransforms as sphcrt
from computeStandardNorms import computeStandardNorms

# %% USER INPUT! ** CHANGE THESE TO RUN **

# RUN: meshes/MakeConvergenceMeshesTempestRemap.sh (Change location of tempestremap executables)
# Change the following variable to the path of the meshes created above
meshesFolder = '/Users/TempestGuerra/Desktop/Remapping-Intercomparison/convergenceMeshes/'

# RUN: ./processConvergenceMeshes.sh
# RUN: ./makeGridConvergenceData_N64_TPO.sh

# %% READ IN THE MESH FILES AND GET COORDINATE DATA

meshRes = ('16', '32', '64', '128')
meshCoords = []
meshCells = []
areas = []
for ff in range(4):
    thisMeshFile = meshesFolder + 'outCSne' + meshRes[ff] + '_enhanced.g'
    thisFid = Dataset(thisMeshFile, 'r')
    meshCoords.append(thisFid.variables['coord'][:])
    meshCells.append(thisFid.variables['connect1'][:])
    areas.append(thisFid.variables['cell_area'][:])
    thisFid.close()
    print('Mesh information for CSne' + meshRes[ff] + '...DONE!')

# %% COMPUTE CENTROID VALUE AND CELL AVERAGES (6th order quadrature) FOR EACH MESH

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
    numCells = len(thisCell)
    NC.append(numCells)

    # Open the sample files for reading
    fnameO1 = 'sample_NM64_O1_outCSne' + str(meshRes[mm]) + '_TPO.nc'
    fidO1 = Dataset(fnameO1, 'r')
    fnameO4 = 'sample_NM64_O8_outCSne' + str(meshRes[mm]) + '_TPO.nc'
    fidO4 = Dataset(fnameO4, 'r')

    # Get the sampled data
    thisValues = fidO1.variables['Topography'][:]
    thisAvgs = fidO4.variables['Topography'][:]

    # Store the differences
    differences.append(np.abs(thisAvgs - thisValues))

    # Store the low order evaluation
    values.append(thisValues)

    # Store the reference (high order cell averages)
    references.append(thisAvgs)

    print('Computed differences for CSne' + meshRes[mm] + '...DONE!')

# %% COMPUTE STANDARD NORMS
ii = 0
L1_test1 = []
L2_test1 = []
Li_test1 = []

for diff in differences:
    # Compute standard norms for tests 1 and 2
    L1_t1, L2_t1, Li_t1 = computeStandardNorms(
        meshCells[ii], values[ii], references[ii], areas[ii], None, False)

    L1_test1.append(L1_t1)
    L2_test1.append(L2_t1)
    Li_test1.append(Li_t1)

    print('L1 Norms by computeStandardNorms()...')
    print(meshRes[ii], L1_test1[ii])
    ii += 1

# %% MAKE THE CONVERGENCE PLOTS (FINGERS CROSSED...)
meshDelta = np.array([1.0, 0.5, 0.25, 0.125])
plt.figure(figsize=(6.0, 10.0))

# L1 Norm plot
scale = L1_test1[0]
order2 = scale * np.power(meshDelta, 2.0)
plt.subplot(3, 1, 1)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, L1_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(
    b=None,
    which='major',
    axis='both',
    color='k',
    linestyle='--',
    linewidth=0.5)
plt.legend(('2nd Order', 'Topography',), loc='upper right')
plt.title('Grid Convergence Test: Centroid Sample to 8th Order Quadrature')
plt.ylabel('L1 Norm Error')

# L2 Norm plot
scale = L2_test1[0]
order2 = scale * np.power(meshDelta, 2.0)
plt.subplot(3, 1, 2)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, L2_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(
    b=None,
    which='major',
    axis='both',
    color='k',
    linestyle='--',
    linewidth=0.5)
plt.legend(('2nd Order', 'Topography'), loc='upper right')
plt.ylabel('L2 Norm Error')

# Li Norm plot
scale = Li_test1[0]
order2 = scale * np.power(meshDelta, 2.0)
plt.subplot(3, 1, 3)
plt.plot(meshDelta, order2, 'k--')
plt.plot(meshDelta, Li_test1, 'bs-')
plt.gca().invert_xaxis()
plt.xscale('log')
plt.yscale('log')
plt.grid(
    b=None,
    which='major',
    axis='both',
    color='k',
    linestyle='--',
    linewidth=0.5)
plt.legend(('2nd Order', 'Topography'), loc='upper right')
plt.xlabel('Normalized Grid Spacing')
plt.ylabel('Linf Norm Error')

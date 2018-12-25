'''
NAME
    NetCDF reader and CANGA intercomparison with Python
PURPOSE
    Reads 3 NetCDF files containing model output (identical variables) and
    computes regridding metrics
PROGRAMMER(S)
    Jorge Guerra, Paul Ullrich
REVISION HISTORY
    
REFERENCES
'''    

import time
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import matplotlib.pyplot as plt
import mpl_toolkits

# Bring in all the different metric modules
from computeGradient import computeGradient
from computeAdjacencyStencil import computeAdjacencyStencil
from computeGlobalConservation import computeGlobalConservation
from computeLocalityMetric import computeLocalityMetric
from computeAreaWeight import computeAreaWeight
from computeStandardNorms import computeStandardNorms
from computeGlobalExtremaMetrics import computeGlobalExtremaMetrics
from computeLocalExtremaMetrics import computeLocalExtremaMetrics

if __name__ == '__main__':
       
       # Set the name of the field variable in question (scalar)
       varName = 'Psi'
       
       # Field sampled at the source (SS)
       nc_fileSS = 'testdata_CSne30_np4_3.nc'
       # Field mapped from source to target (S2T)
       nc_fileS2T = 'testdata_CSne30_2_RLL1deg_np4_3.nc'
       # Field sampled at the target (ST)
       nc_fileST = 'testdata_RLL1deg_np4_3.nc'
       
       # Source Exodus .g file (only needed for global conservation)
       exo_fileS = 'outCSne30.g'
       # Target Exodus .g file
       exo_fileT = 'outRLL1deg.g'
       
       # Open the .g mesh files for reading
       g_fidS = Dataset(exo_fileS)
       g_fidT = Dataset(exo_fileT)
       
       # Get connectivity and coordinate arrays
       varConS = g_fidS.variables['connect1'][:]
       varCoordS = g_fidS.variables['coord'][:]
       varConT = g_fidT.variables['connect1'][:]
       varCoordT = g_fidT.variables['coord'][:]
       
       start = time.time()
       # Compute adjacency maps for both meshes (source stencil NOT needed)
       #edgeMapS, sortedEdgeMapS, varConStenDexS = computeAdjacencyStencil(varConS) 
       edgeMapT, sortedEdgeMapT, varConStenDexT = computeAdjacencyStencil(varConT)
       endt = time.time()
       print('Time to precompute adjacency maps (sec): ', endt - start)
       
       start = time.time()
       # Precompute the area weights and then look them up in the integral below
       NEL = len(varConS)
       areaS = np.zeros((NEL,1))
       for ii in range(NEL):
              areaS[ii] = computeAreaWeight(varCoordS, varConS[ii, :])
              
       # Precompute the area weights and then look them up in the integrals below
       NEL = len(varConT)
       areaT = np.zeros((NEL,1))
       for ii in range(NEL):
              areaT[ii] = computeAreaWeight(varCoordT, varConT[ii, :])
       
       endt = time.time()
       print('Time to precompute mesh areas (sec): ', endt - start)
       
       start = time.time()
       # Open the .nc data files for reading
       nc_fidSS = Dataset(nc_fileSS, 'r')
       nc_fidS2T = Dataset(nc_fileS2T, 'r')
       nc_fidST = Dataset(nc_fileST, 'r')
       
       # Get the SS data
       varSS = nc_fidSS.variables[varName][:]
       # Get the S2T data
       varS2T = nc_fidS2T.variables[varName][:]
       # Get the ST data
       varST = nc_fidST.variables[varName][:]
       
       # Check the extracted variables for dimensionality
       # If the variables are 2D then reshape along the longitude (ASSUMED)
       # TO DO: Check with Paul and tempestremap code on reshaping
       VS = varSS.shape
       if len(VS) > 1:
              varSS = np.reshape(varSS, VS[0] * VS[1])
              
       VS = varS2T.shape
       if len(VS) > 1:
              varS2T = np.reshape(varS2T, VS[0] * VS[1])
              
       VS = varST.shape
       if len(VS) > 1:
              varST = np.reshape(varST, VS[0] * VS[1])
       
       endt = time.time()
       print('Time to read NC and Exodus data (sec): ', endt - start)
       
       start = time.time()
       # Precompute the gradient operator on regridded and sampled target data
       varsOnTargetMesh = [varST, varS2T]
       gradientST = computeGradient(varsOnTargetMesh, varCoordT, varConStenDexT, areaT)
       
       endt = time.time()
       print('Time to compute gradients on target mesh (sec): ')
       
       start = time.time()
       # Global conservation metric
       L_g = computeGlobalConservation(varSS, varS2T, varST, areaS, areaT)
       # Locality measure (returns an array for each target DOF)
       #L_local = computeLocalityMetric(varS2T, varST, varConT, varCoordT)
       # Standard Error norms (L_1, L_2, L_inf)
       L_1, L_2, L_inf = computeStandardNorms(varS2T, varST, areaT)
       # Global Extrema preservation
       Lmin, Lmax = computeGlobalExtremaMetrics(varS2T, varST)
       # Local Extrema preservation
       Lmin_1, Lmin_2, Lmin_inf, Lmax_1, Lmax_2, Lmax_inf = \
       computeLocalExtremaMetrics(areaT, varSS, varS2T, varST, varConS, varCoordS, varConT, varCoordT)
       endt = time.time()
       print('Time to execute metrics (sec): ', endt - start)
       
       # Make some plots here
       
       # Close original NetCDF file.
       nc_fidSS.close()
       nc_fidS2T.close()
       nc_fidST.close()

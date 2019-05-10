#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
NAME
    NetCDF reader and CANGA intercomparison with Python
PURPOSE
    Reads 3 NetCDF files containing model output (identical variables) and
    computes regridding metrics. Also takes mesh data from Exodus or SCRIP.
PROGRAMMER(S)
    Jorge Guerra, Paul Ullrich
REVISION HISTORY
    
REFERENCES
'''    
#%%
import sys, getopt
import time
import math as mt
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

# Bring in all the different metric modules
from computeGradient2 import computeGradient2
from computeGlobalConservation import computeGlobalConservation
#from computeLocalityMetric import computeLocalityMetric
from computeStandardNorms import computeStandardNorms
from computeGlobalExtremaMetrics import computeGlobalExtremaMetrics
from computeLocalExtremaMetrics import computeLocalExtremaMetrics
from computeGradientPreserveMetrics import computeGradientPreserveMetrics

def computeLL2Cart(cellCoord):
       # Loop over the Lon/Lat coordinate array, extract Cartesian coords
       # Input array is [lon, lat, radius] (radians)
       NC = np.size(cellCoord, axis=0)
       varCart = np.zeros((NC, 3))
       for ii in range(NC):
              RO = cellCoord[ii,2]
              lon = cellCoord[ii,0]
              lat = cellCoord[ii,1]
              X = RO * mt.cos(lat) * mt.sin(lon)
              Y = RO * mt.cos(lat) * mt.cos(lon)
              Z = RO * mt.sin(lat)
              RC = mt.sqrt(X**2 + Y**2 + Z**2)
              varCart[ii,:] = [X, Y, Z]
              varCart[ii,:] *= 1.0 / RC
       
       # INPUT IS IN RADIANS
       return varCart

# Parse the command line
def parseCommandLine(argv):
       # Field variable name
       varName = ''
       
       # NC files with field data
       sourceSampledFile = ''
       targetSampledFile = ''
       remappedFile = ''
       
       # Mesh information files
       sourceMesh = ''
       targetMesh = ''
       
       ExodusSingleConn = False
       SCRIPwithoutConn = False
       SCRIPwithConn = False
       
       # Mesh data configuration
       meshConfig = ''
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['ss=','s2t=','st=','sm=','tm=', \
                                         'ExodusSingleConn', 'SCRIPwithoutConn', \
                                         'SCRIPwithConn'])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAMEtricsDriver.py', \
                    '-ss <sourceSampledFile>', \
                    '-s2t <remappedFile>', \
                    '-st <targetSampledFile>', \
                    '-sm <sourceMesh>', \
                    '-tm <targetMesh>', \
                    '--<meshConfiguration>')
              sys.exit(2)
              
       for opt, arg in opts:
              # Request for usage help
              if opt == '-h':
                     print('Command line not properly set:', \
                           'CANGAMEtricsDriver.py', \
                           '-ss <sourceSampledFile>', \
                           '-s2t <remappedFile>', \
                           '-st <targetSampledFile>', \
                           '-sm <sourceMesh>', \
                           '-tm <targetMesh>', \
                           '--<meshConfiguration>')
                     sys.exit()
              elif opt == '-v':
                     varName = arg
              elif opt == '--ss':
                     sourceSampledFile = arg
              elif opt == '--s2t':
                     remappedFile = arg
              elif opt == '--st':
                     targetSampledFile = arg
              elif opt == '--sm':
                     sourceMesh = arg
              elif opt == '--tm':
                     targetMesh = arg
              elif opt == '--ExodusSingleConn':
                     meshConfig = 'ExodusSingleConn'
                     ExodusSingleConn = True
              elif opt == '--SCRIPwithoutConn':
                     meshConfig = 'SCRIPwithoutConn'
                     SCRIPwithoutConn = True
              elif opt == '--SCRIPwithConn':
                     meshConfig = 'SCRIPwithConn'
                     SCRIPwithConn = True
                     
       # Check that mesh files are consistent
       if sourceMesh[len(sourceMesh) - 1] != targetMesh[len(targetMesh) - 1]:
              print('Mesh data files have inconsistent file extensions!')
              print('Source and target mesh files MUST have the same extension.')
              sys.exit(2)
              
       # Check that only one configuration is chosen
       if (ExodusSingleConn == True) & (SCRIPwithoutConn == True):
              print('Expecting only ONE mesh configuration option!')
              print('Multiple options are set.')
              sys.exit(2)
       elif (ExodusSingleConn == True) & (SCRIPwithConn == True):
              print('Expecting only ONE mesh configuration option!')
              print('Multiple options are set.')
              sys.exit(2)
       elif (SCRIPwithoutConn == True) & (SCRIPwithConn == True):
              print('Expecting only ONE mesh configuration option!')
              print('Multiple options are set.')
              sys.exit(2)
       elif (ExodusSingleConn == True) & (SCRIPwithoutConn == True) & (SCRIPwithConn == True):
              print('Expecting only ONE mesh configuration option!')
              print('Multiple options are set.')
              sys.exit(2)
       elif (ExodusSingleConn == False) & (SCRIPwithoutConn == False) & (SCRIPwithConn == False):
              print('ONE mesh configuration option must be set!')
              print('None of the options are set.')
              sys.exit(2)
                     
       # Check that configurations match the correct file extensions
       if meshConfig == 'ExodusSingleConn':
              if sourceMesh[len(sourceMesh) - 1] != 'g':
                     print('Expected Exodus .g mesh data files!')
                     print('Exodus files MUST have .g extension')
                     sys.exit(2)
       if meshConfig == 'SCRIPwithoutConn':
              if sourceMesh[len(sourceMesh) - 1] != 'c':
                     print('Expected SCRIP .nc mesh data files!')
                     print('SCRIP files MUST have .nc extension')
                     sys.exit(2)
       
       return varName, \
              sourceSampledFile, remappedFile, targetSampledFile, \
              sourceMesh, targetMesh, \
              ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison metrics!')
       print('Authors: Jorge Guerra, Paul Ullrich, 2019')
       # Global parameters
       #kdleafs = 100
       #COINCIDENT_TOLERANCE = 1.0E-14

       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       varName, nc_fileSS, nc_fileS2T, nc_fileST, mesh_fileS, mesh_fileT, \
       ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn = \
       parseCommandLine(sys.argv[1:])
       
       # Set the names for the auxiliary area and adjacency maps (NOT USER)
       varAreaName = 'cell_area'
       varAdjaName = 'cell_edge_adjacency'
       
       """ SET INPUT HERE FOR DEVELOPMENT TESTING
       # Set the mesh configuration (mutually exclusive):
       # ExodusSingleConn -> DEFAULT BEST (DEGENERATE POLYGONS OK)
       # ExodusMultiConn -> NOT IMPLEMENTED (DEGENERATE POLYGONS OK)
       # SCRIPwithoutConn -> UNFORTUNATE SECOND
       # SCRIPwithConn -> NOT IMPLEMENTED (READING METIS MESH INFO PROBLEMATIC)
       ExodusSingleConn = True
       #ExodusMultiConn = False
       SCRIPwithoutConn = False
       SCRIPwithConn = False
       
       if ExodusSingleConn:
              # Source Exodus .g file
              mesh_fileS = 'outCSne30.g'
              # Target Exodus .g file
              mesh_fileT = 'outRLL1deg.g'
              #mesh_fileT = 'outICO64.g'
              
       if SCRIPwithoutConn:
              # Source SCRIP file
              mesh_fileS = 'Grids/ne30np4_pentagons.091226.nc'
              # Target SCRIP file
              mesh_fileT = 'Grids/ne30np4_latlon.091226.nc'
       
       # Set the name of the field variable in question (scalar)
       #varName = 'Psi'
       varName = 'TotalPrecipWater'
       #varName = 'CloudFraction'
       #varName = 'Topography'
       
       # Field sampled at the source (SS)
       nc_fileSS = 'testdata_outCSne30_TPW_CFR_TPO.nc'
       #nc_fileSS = 'testdata_CSne30_np4_3.nc'
       # Field mapped from source to target (S2T)
       #nc_fileS2T = 'testdata_CSne30_2_RLL1deg_np4_3.nc'
       #nc_fileS2T = 'testdata_CSne30_2_ICO64_np4_3.nc'
       nc_fileS2T = 'testdata_outCSne30_2_RLL1deg_TPW.nc'
       #nc_fileS2T = 'testdata_outCSne30_2_RLL1deg_CFR.nc'
       #nc_fileS2T = 'testdata_outCSne30_2_RLL1deg_TPO.nc'
       #nc_fileS2T = 'testdata_outCSne30_2_ICO64_TPW.nc'
       #nc_fileS2T = 'testdata_outCSne30_2_ICO64_CFR.nc'
       #nc_fileS2T = 'testdata_outCSne30_2_ICO64_TPO.nc'
       # Field sampled at the target (ST)
       nc_fileST = 'testdata_outRLL1deg_TPW_CFR_TPO.nc'
       #nc_fileST = 'testdata_outICO64_TPW_CFR_TPO.nc'
       #nc_fileST = 'testdata_RLL1deg_np4_3.nc'
       #nc_fileST = 'testdata_ICO64_np4_3.nc'
       #""" 
       
       if ExodusSingleConn:
              numEdges = 'num_nod_per_el1'
              numCells = 'num_el_in_blk1'
              numDims = 'cart_dims'
              
              # Open the .g mesh files for reading
              m_fidS = Dataset(mesh_fileS, 'r')
              m_fidT = Dataset(mesh_fileT, 'r')
              
              # Get connectivity and coordinate arrays (check for multiple connectivity)
              varConS = m_fidS.variables['connect1'][:]
              varCoordS = m_fidS.variables['coord'][:]
              varConT = m_fidT.variables['connect1'][:]
              varCoordT = m_fidT.variables['coord'][:]
              
       elif SCRIPwithoutConn:
              numEdges = 'grid_corners'
              numCells = 'grid_size'
              numDims = 'cart_dims'
              connCell = 'element_corners_id'
              coordCell = 'grid_corners_cart'
              numVerts = 'grid_corners_size'
              
              # Open the .nc SCRIP files for reading
              m_fidS = Dataset(mesh_fileS, 'a')
              m_fidT = Dataset(mesh_fileT, 'a')
                    
              start = time.time()
              try:
                     print('Reading coordinate and connectivity from augmented SCRIP')
                     varConS = m_fidS.variables[connCell][:]
                     varConT = m_fidT.variables[connCell][:]
                     varCoordS = m_fidS.variables[coordCell][:]
                     varCoordT = m_fidT.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS MESH FILE!')
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
       elif SCRIPwithConn:
              numEdges = 'ncorners'
              numCells = 'ncells'
              numDims = 'cart_dims'
              connCell = 'element_corners'
              coordCell = 'grid_corners_cart'
              
              # Open the .nc SCRIP files for reading
              m_fidS = Dataset(mesh_fileS, 'a')
              m_fidT = Dataset(mesh_fileT, 'a')
              
              # Get connectivity and coordinate arrays
              varConS = m_fidS.variables[connCell][:]
              varConS = varConS.T
              varConT = m_fidT.variables[connCell][:]
              varConT = varConT.T
                     
              start = time.time()
              try:
                     print('Reading coordinate and connectivity from augmented SCRIP')
                     varCoordS = m_fidS.variables[coordCell][:]
                     varCoordT = m_fidT.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS MESH FILE!')
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
              
       m_fidS.close()
       m_fidT.close()
       
       
       #%%
       start = time.time()
       print('Reading adjacency maps...')
       # Store the adjacency map in the original grid netcdf file (target mesh)
       m_fidT = Dataset(mesh_fileT, 'a')
       # Check for existing variable data
       try:
              if m_fidT.variables[varAdjaName].name == varAdjaName:
                     varConStenDexT = m_fidT.variables[varAdjaName][:]
                     
       except:
              print('PRE-PROCESSING FOR ADJACENCY NOT DONE ON THIS MESH FILE!')
       m_fidT.close()
              
       endt = time.time()
       print('Time to read adjacency maps (sec): ', endt - start)
       #%%
       start = time.time()
       print('Reading source and target mesh areas...')
       
       # Store the grid cell areas in the original netcdf file (source and target)
       m_fidS = Dataset(mesh_fileS, 'a')
       # Check for existing variable data
       try:
              if m_fidS.variables[varAreaName].name == varAreaName:
                     areaS = m_fidS.variables[varAreaName][:]
       
       except:
              print('PRE-PROCESSING FOR AREAS NOT DONE ON SOURCE MESH FILE!')
              
       m_fidS.close()
              
       m_fidT = Dataset(mesh_fileT, 'a')
       # Check for existing variable data
       try:
              if m_fidT.variables[varAreaName].name == varAreaName:
                     areaT = m_fidT.variables[varAreaName][:]
       
       except:
              print('PRE-PROCESSING FOR AREAS NOT DONE ON TARGET MESH FILE!')
              
       m_fidT.close()
       
       endt = time.time()
       print('Time to precompute/read mesh areas (sec): ', endt - start)

       #%%
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
       VS = varSS.shape
       if len(VS) > 1:
              varSS = np.reshape(varSS, VS[0] * VS[1])
              
       VS = varS2T.shape
       if len(VS) > 1:
              varS2T = np.reshape(varS2T, VS[0] * VS[1])
              
       VS = varST.shape
       if len(VS) > 1:
              varST = np.reshape(varST, VS[0] * VS[1])
              
       #%% Close original NetCDF file.
       nc_fidSS.close()
       nc_fidS2T.close()
       nc_fidST.close()
       
       endt = time.time()
       print('Time to read NC and Exodus data (sec): ', endt - start)
       #%%
       start = time.time()
       print('Computing or reading gradients for target sampled and regridded fields...')
       
       # Open data files for storage of gradient data
       nc_fidS2T = Dataset(nc_fileS2T, 'a')
       nc_fidST = Dataset(nc_fileST, 'a')
       
       varGradientName = 'FieldGradient'
       try:
              # Read in previously stored data if it exists
              if nc_fidST.variables[varGradientName].name == varGradientName:
                     gradST = nc_fidST.variables[varGradientName][:]
                     
              if nc_fidS2T.variables[varGradientName].name == varGradientName:
                     gradS2T = nc_fidS2T.variables[varGradientName][:]
                     
              varsOnTM = [varST, varS2T]
              gradientsOnTM = [gradST, gradS2T]
       except KeyError:
              # Precompute the gradients on target mesh ONLY once
              varsOnTM = [varST, varS2T]
              gradientsOnTM, cellCoordT = computeGradient2(varsOnTM, varConT, varCoordT, varConStenDexT, areaT)
              
              # Create new Cartesian Earth centered vector dimensions
              try:
                     nc_fidST.createDimension(numDims, 3)
              except RuntimeError:
                     print('Dimensions for gradient variable already exist in field data file.')
              try:
                     nc_fidS2T.createDimension(numDims, 3)
              except RuntimeError:
                     print('Dimensions for gradient variable already exist in field data file.')
              
              # Create new dimension for the number of cells
              try:
                     nc_fidST.createDimension(numCells, np.size(varST, axis=0))
              except RuntimeError:
                     print('Dimensions for gradient variable already exist in field data file.')
              try:
                     nc_fidS2T.createDimension(numCells, np.size(varS2T, axis=0))
              except RuntimeError:
                     print('Dimensions for gradient variable already exist in field data file.')
              
              # Store the gradients on target mesh
              try:
                     gradFileOut = nc_fidST.createVariable(varGradientName, 'f8', (numDims, numCells))
                     gradFileOut[:] = gradientsOnTM[0]
              except RuntimeError:
                     print('Gradient variable already exists in ST field data file.')
              
              try:
                     gradFileOut = nc_fidS2T.createVariable(varGradientName, 'f8', (numDims, numCells))
                     gradFileOut[:] = gradientsOnTM[1]
              except RuntimeError:
                     print('Gradient variable already exists in S2T field data file.')
                     
       
       nc_fidS2T.close()
       nc_fidST.close()
       endt = time.time()
       print('Time to compute/read gradients on target mesh (sec): ', endt - start)
       
       #%%
       start = time.time()
       print('Computing all metrics...')
       # Global conservation metric
       massSS, massS2T, massST, L_g = computeGlobalConservation(varSS, varS2T, varST, areaS, areaT)
       # Locality measure (returns an array for each target DOF)
       #L_local = computeLocalityMetric(varS2T, varST, varConT, varCoordT)
       # Standard Error norms (L_1, L_2, L_inf)
       L_1, L_2, L_inf = computeStandardNorms(varS2T, varST, areaT)
       # Global Extrema preservation
       Lmin, Lmax = computeGlobalExtremaMetrics(varS2T, varST)
       # Local Extrema preservation
       Lmin_1, Lmin_2, Lmin_inf, Lmax_1, Lmax_2, Lmax_inf = \
       computeLocalExtremaMetrics(areaT, varSS, varS2T, varST, varConS, varCoordS, varConT, varCoordT)
       # Gradient preservation
       H1, H1_2 = computeGradientPreserveMetrics(gradientsOnTM, varsOnTM, areaT)
       endt = time.time()
       print('Time to execute metrics (sec): ', endt - start)
       #%%
       # Print out a table with metric results
       print('Global conservation: %16.15e' % np.ravel(L_g))
       print('Global L1 error:     %16.15e' % np.ravel(L_1))
       print('Global L2 error:     %16.15e' % np.ravel(L_2))
       print('Global Linf error:   %16.15e' % np.ravel(L_inf))
       print('Global max error:    %16.15e' % np.ravel(Lmax))
       print('Global min error:    %16.15e' % np.ravel(Lmin))
       print('Local max L1 error:  %16.15e' % np.ravel(Lmax_1))
       print('Local max L2 error:  %16.15e' % np.ravel(Lmax_2))
       print('Local max Lm error:  %16.15e' % np.ravel(Lmax_inf))
       print('Local min L1 error:  %16.15e' % np.ravel(Lmin_1))
       print('Local min L2 error:  %16.15e' % np.ravel(Lmin_2))
       print('Local min Lm error:  %16.15e' % np.ravel(Lmin_inf))
       print('Gradient semi-norm:  %16.15e' % np.ravel(H1_2))
       print('Gradient full-norm:  %16.15e' % np.ravel(H1))       
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
from computeGradientSE import computeGradientSE
from computeGradientFV2 import computeGradientFV2
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
       targetSampledFile = ''
       remappedFile = ''
       
       # Mesh information files
       targetMesh = ''
       
       ExodusSingleConn = False
       SCRIPwithoutConn = False
       SCRIPwithConn = False
       SpectralElement = False
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['s2t=', 'st=', 'tm=', \
                                         'ExodusSingleConn', 'SCRIPwithoutConn', \
                                         'SCRIPwithConn', 'SpectralElement'])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAMEtricsDriver.py', \
                    '-s2t <remappedFile>', \
                    '-st <targetSampledFile>', \
                    '-tm <targetMesh>', \
                    '--<meshConfiguration>', \
                    '--<isSpectralElementMesh>')
              sys.exit(2)
              
       for opt, arg in opts:
              # Request for usage help
              if opt == '-h':
                     print('Command line not properly set:', \
                           'CANGAMEtricsDriver.py', \
                           '-s2t <remappedFile>', \
                           '-st <targetSampledFile>', \
                           '-tm <targetMesh>', \
                           '--<meshConfiguration>', \
                           '--<isSpectralElementMesh>')
                     sys.exit()
              elif opt == '-v':
                     varName = arg
              elif opt == '--s2t':
                     remappedFile = arg
              elif opt == '--st':
                     targetSampledFile = arg
              elif opt == '--tm':
                     targetMesh = arg
              elif opt == '--ExodusSingleConn':
                     ExodusSingleConn = True
              elif opt == '--SCRIPwithoutConn':
                     SCRIPwithoutConn = True
              elif opt == '--SCRIPwithConn':
                     SCRIPwithConn = True
              elif opt == '--SpectralElement':
                     SpectralElement = True
              
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
       
       return varName, remappedFile, targetSampledFile, targetMesh, \
              ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn, SpectralElement

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison metrics!')
       print('Authors: Jorge Guerra, Paul Ullrich, 2019')

       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       varName, nc_fileS2T, nc_fileST, mesh_fileT, \
       ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn, SpectralElement = \
       parseCommandLine(sys.argv[1:])
       
       # Set the names for the auxiliary area and adjacency maps (NOT USER)
       varAreaName = 'cell_area'
       varJacoName = 'element_jacobians'
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
              
              if SpectralElement:
                     connCell = 'element_gll_conn'
                     coordCell = 'grid_gll_cart'
              else:
                     connCell = 'connect1'
                     coordCell = 'coord'
              
              # Open the .g mesh files for reading
              m_fidT = Dataset(mesh_fileT, 'r')
              
              # Get connectivity and coordinate arrays (check for multiple connectivity)
              varConT = m_fidT.variables[connCell][:]
              varCoordT = m_fidT.variables[coordCell][:]
              
       elif SCRIPwithoutConn:
              numEdges = 'grid_corners'
              numCells = 'grid_size'
              numDims = 'cart_dims'
              numVerts = 'grid_corners_size'
              
              if SpectralElement:
                     connCell = 'element_gll_conn'
                     coordCell = 'grid_gll_cart'
              else:
                     connCell = 'element_corners_id'
                     coordCell = 'grid_corners_cart'
              
              # Open the .nc SCRIP files for reading
              m_fidT = Dataset(mesh_fileT, 'a')
                    
              start = time.time()
              try:
                     print('Reading coordinate and connectivity from augmented SCRIP')
                     varConT = m_fidT.variables[connCell][:]
                     varCoordT = m_fidT.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS MESH FILE!')
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
       elif SCRIPwithConn:
              numEdges = 'ncorners'
              numCells = 'ncells'
              numDims = 'cart_dims'
              
              if SpectralElement:
                     connCell = 'element_gll_conn'
                     coordCell = 'grid_gll_cart'
              else:
                     connCell = 'element_corners'
                     coordCell = 'grid_corners_cart' 
              
              # Open the .nc SCRIP files for reading
              m_fidT = Dataset(mesh_fileT, 'a')
              
              # Get connectivity and coordinate arrays
              varConT = m_fidT.variables[connCell][:]
              varConT = varConT.T
                     
              start = time.time()
              try:
                     print('Reading coordinate and connectivity from augmented SCRIP')
                     varCoordT = m_fidT.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS MESH FILE!')
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
              
       m_fidT.close()
       
       
       #%%
       start = time.time()
       print('Reading adjacency maps...')
       # Fetch the adjacency map in the original grid netcdf file (target mesh)
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
       print('Reading mesh areas...')
              
       m_fidT = Dataset(mesh_fileT, 'a')
       # Check for existing variable data
       try:
              if m_fidT.variables[varAreaName].name == varAreaName:
                     areaT = m_fidT.variables[varAreaName][:]
       
       except:
              print('PRE-PROCESSING FOR AREAS NOT DONE ON TARGET MESH FILE!')
              
       m_fidT.close()
       
       endt = time.time()
       print('Time to read mesh areas (sec): ', endt - start)
       
       #%%
       if SpectralElement:
              start = time.time()
              print('Reading SE mesh jacobians...')
                     
              m_fidT = Dataset(mesh_fileT, 'a')
              # Check for existing variable data
              try:
                     if m_fidT.variables[varJacoName].name == varJacoName:
                            jacobiansT = m_fidT.variables[varJacoName][:]
              
              except:
                     print('PRE-PROCESSING FOR JACOBIANS NOT AVAILABLE ON TARGET MESH FILE!')
                     
              m_fidT.close()
              
              endt = time.time()
              print('Time to read SE mesh jacobians (sec): ', endt - start)
       else:
              jacobiansT = None

       #%%
       start = time.time()
       # Open the .nc data files for reading
       nc_fidS2T = Dataset(nc_fileS2T, 'r')
       nc_fidST = Dataset(nc_fileST, 'r')
       
       # Get the S2T data
       varS2T = nc_fidS2T.variables[varName][:]
       # Get the ST data
       varST = nc_fidST.variables[varName][:]
       
       # Check the extracted variables for dimensionality
       # If the variables are 2D then reshape along the longitude (ASSUMED)     
       VS = varS2T.shape
       if len(VS) > 1:
              varS2T = np.reshape(varS2T, VS[0] * VS[1])
              
       VS = varST.shape
       if len(VS) > 1:
              varST = np.reshape(varST, VS[0] * VS[1])
              
       #%% Close original NetCDF file.
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
              #gradientsOnTM, cellCoordT = computeGradient2(varsOnTM, varConT, varCoordT, varConStenDexT, areaT)
       except KeyError:
              # Precompute the gradients on target mesh ONLY once
              varsOnTM = [varST, varS2T]
              
              if SpectralElement:
                     numDOFS = coordCell
                     gradientsOnTM = computeGradientSE(varsOnTM, varConT, varCoordT, jacobiansT)
              else: 
                     numDOFS = numCells
                     gradientsOnTM = computeGradientFV2(varsOnTM, varConT, varCoordT, varConStenDexT, areaT)
              
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
                     nc_fidST.createDimension(numDOFS, np.size(varST, axis=0))
              except RuntimeError:
                     print('Dimensions for gradient variable already exist in field data file.')
              try:
                     nc_fidS2T.createDimension(numDOFS, np.size(varS2T, axis=0))
              except RuntimeError:
                     print('Dimensions for gradient variable already exist in field data file.')
              
              # Store the gradients on target mesh
              try:
                     gradFileOut = nc_fidST.createVariable(varGradientName, 'f8', (numDims, numDOFS))
                     gradFileOut[:] = gradientsOnTM[0]
              except RuntimeError:
                     print('Gradient variable already exists in ST field data file.')
              
              try:
                     gradFileOut = nc_fidS2T.createVariable(varGradientName, 'f8', (numDims, numDOFS))
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
       massS2T, massST, L_g = computeGlobalConservation(varConT, varS2T, varST, areaT, jacobiansT, SpectralElement)
       # Locality measure (returns an array for each target DOF)
       #L_local = computeLocalityMetric(varS2T, varST, varConT, varCoordT)
       # Standard Error norms (L_1, L_2, L_inf)
       L_1, L_2, L_inf = computeStandardNorms(varConT, varS2T, varST, areaT, jacobiansT, SpectralElement)
       # Global Extrema preservation
       Lmin, Lmax = computeGlobalExtremaMetrics(varS2T, varST)
       # Local Extrema preservation
       Lmin_inf, Lmax_inf = \
       computeLocalExtremaMetrics(varConStenDexT, varConT, varCoordT, varS2T, varST, SpectralElement)
       # Gradient preservation
       H1, H1_2 = computeGradientPreserveMetrics(varConT, gradientsOnTM, varsOnTM, areaT, jacobiansT, SpectralElement)
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
       print('Local max Lm error:  %16.15e' % np.ravel(Lmax_inf))
       print('Local min Lm error:  %16.15e' % np.ravel(Lmin_inf))
       print('Gradient semi-norm:  %16.15e' % np.ravel(H1_2))
       print('Gradient full-norm:  %16.15e' % np.ravel(H1))       
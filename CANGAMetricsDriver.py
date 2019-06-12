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
from computeGradientFV3 import computeGradientFV3
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
       sourceMeshFile = ''
       targetMeshFile = ''
       sourceMeshConfig = 0
       targetMeshConfig = 0

       sourceSE = False
       targetSE = False
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['ss=', 's2t=', 'st=', 'sm=', 'tm=', \
                                         'smc=', 'tmc=', 'sourceSE', 'targetSE'])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAMEtricsDriver.py', \
                    '--ss <SourceSampledFile>', \
                    '--s2t <remappedFile>', \
                    '--st <targetSampledFile>', \
                    '--sm <sourceMeshFile>', \
                    '--smc <sourceMeshConfiguration>', \
                    '--<isSourceSpectralElementMesh>', \
                    '--tm <targetMeshFile>', \
                    '--tmc <targetMeshConfiguration>', \
                    '--<isTargetSpectralElementMesh>')
              sys.exit(2)
              
       for opt, arg in opts:
              # Request for usage help
              if opt == '-h':
                     print('Command line not properly set:', \
                           'CANGAMEtricsDriver.py', \
                           '--ss <SourceSampledFile>', \
                           '--s2t <remappedFile>', \
                           '--st <targetSampledFile>', \
                           '--sm <sourceMeshFile>', \
                           '--smc <sourceMeshConfiguration>', \
                           '--<isSourceSpectralElementMesh>', \
                           '--tm <targetMeshFile>', \
                           '--tmc <targetMeshConfiguration>', \
                           '--<isTargetSpectralElementMesh>')
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
                     sourceMeshFile = arg
              elif opt == '--tm':
                     targetMeshFile = arg
              elif opt == '--smc':
                     sourceMeshConfig = int(arg)
              elif opt == '--tmc':
                     targetMeshConfig = int(arg)
              elif opt == '--sourceSE':
                     sourceSE = True
              elif opt == '--targetSE':
                     targetSE = True
                                   
       # Input checks
       if sourceMeshConfig > 3:
              print('ERROR: Invalid source mesh configuration (1-3)')
              sys.exit(2)
       
       if targetMeshConfig > 3:
              print('ERROR: Invalid target mesh configuration (1-3)')
              sys.exit(2)
       
       return varName, sourceSampledFile, remappedFile, targetSampledFile, \
              sourceMeshFile, targetMeshFile, \
              sourceMeshConfig, targetMeshConfig, sourceSE, targetSE
              
def loadMeshData(mesh_file, mesh_config, SpectralElement):
       
       if mesh_config == 1:
              numEdges = 'num_nod_per_el1'
              numCells = 'num_el_in_blk1'
              numDims = 'cart_dims'
              numVerts = ''
              
              if SpectralElement:
                     connCell = 'element_gll_conn'
                     coordCell = 'grid_gll_cart'
              else:
                     connCell = 'connect1'
                     coordCell = 'coord'
              
              # Open the .g mesh files for reading
              m_fid = Dataset(mesh_file, 'a')
              
              # Get connectivity and coordinate arrays (check for multiple connectivity)
              varConn = m_fid.variables[connCell][:]
              varCoord = m_fid.variables[coordCell][:]
              
       elif mesh_config == 2:
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
              m_fid = Dataset(mesh_file, 'a')
                    
              start = time.time()
              try:
                     print('Reading coordinate and connectivity from augmented SCRIP')
                     varConn = m_fid.variables[connCell][:]
                     varCoord = m_fid.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS SCRIP MESH FILE!')
                     sys.exit()
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
       
       elif mesh_config == 3:
              numEdges = 'ncorners'
              numCells = 'ncells'
              numDims = 'cart_dims'
              numVerts = ''
              
              if SpectralElement:
                     connCell = 'element_gll_conn'
                     coordCell = 'grid_gll_cart'
              else:
                     connCell = 'element_corners'
                     coordCell = 'grid_corners_cart' 
              
              # Open the .nc SCRIP files for reading
              m_fid = Dataset(mesh_file, 'a')
              
              # Get connectivity and coordinate arrays
              varConn = m_fid.variables[connCell][:]
              varConn = varConn.T
                     
              start = time.time()
              try:
                     print('Reading coordinate and connectivity from augmented SCRIP')
                     varCoord = m_fid.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS SCRIP MESH FILE!')
                     sys.exit()
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
              
       m_fid.close()
       
       return varCoord, varConn, numEdges, numCells, numDims, numVerts

def loadMeshAreas(mesh_file, varAreaName):
       start = time.time()
       print('Reading mesh areas...')
              
       m_fid = Dataset(mesh_file, 'a')
       # Check for existing variable data
       try:
              if m_fid.variables[varAreaName].name == varAreaName:
                     areas = m_fid.variables[varAreaName][:]
       
       except:
              print('PRE-PROCESSING FOR AREAS NOT DONE ON TARGET MESH FILE!')
              sys.exit()
              
       m_fid.close()
       
       endt = time.time()
       print('Time to read mesh areas (sec): ', endt - start)
       
       return areas

def loadMeshAdjacencyMap(mesh_file, varAdjaName):
       start = time.time()
       print('Reading adjacency maps...')
       # Fetch the adjacency map in the original grid netcdf file (target mesh)
       m_fid = Dataset(mesh_file, 'a')
       # Check for existing variable data
       try:
              if m_fid.variables[varAdjaName].name == varAdjaName:
                     varConStenDex = m_fid.variables[varAdjaName][:]
                     
       except:
              print('PRE-PROCESSING FOR ADJACENCY NOT DONE ON THIS MESH FILE!')
              sys.exit()
              
       m_fid.close()
              
       endt = time.time()
       print('Time to read adjacency maps (sec): ', endt - start)
       
       return varConStenDex

def loadMeshJacobians(mesh_file, varJacoName, SpectralElement):
       if SpectralElement:
              start = time.time()
              print('Reading SE mesh jacobians...')
                     
              m_fid = Dataset(mesh_file, 'a')
              # Check for existing variable data
              try:
                     if m_fid.variables[varJacoName].name == varJacoName:
                            jacobians = m_fid.variables[varJacoName][:]
              
              except:
                     print('ERROR: PRE-PROCESSING FOR JACOBIANS NOT AVAILABLE ON TARGET MESH FILE!')
                     sys.exit()
                     
              m_fid.close()
              
              endt = time.time()
              print('Time to read SE mesh jacobians (sec): ', endt - start)
       else:
              jacobians = None
              
       return jacobians

def loadField(var_file, varName):
       start = time.time()
       # Open the .nc data files for reading
       nc_fid = Dataset(var_file, 'a')
       
       # Get the field data
       varField = nc_fid.variables[varName][:]
       
       # Check the extracted variables for dimensionality
       # If the variables are 2D then reshape along the longitude (ASSUMED)     
       VS = varField.shape
       if len(VS) > 1:
              varField = np.reshape(varField, VS[0] * VS[1])
              
       #%% Close original NetCDF file.
       nc_fid.close()
       
       endt = time.time()
       print('Time to read NC and Exodus data (sec): ', endt - start)
       
       return varField

def loadFieldGradient(var_file, varGradientName, varField, varConn, varCoord, varConStenDex, jacobians, numCells, numDims, SpectralElement):
       
       start = time.time()
       print('Computing or reading gradients for target sampled and regridded fields...')
       
       # Open data files for storage of gradient data
       nc_fid = Dataset(var_file, 'a')
        
       # Read in previously stored ST data if it exists, or compute it and store
       try:
              if nc_fid.variables[varGradientName].name == varGradientName:
                     gradField = nc_fid.variables[varGradientName][:]
       except KeyError:
              if SpectralElement:
                     # This comes from mesh preprocessing
                     numDOFS = 'grid_gll_size'
                     gradField = computeGradientSE(varField, varConn, varCoord, 4, jacobians)
              else: 
                     numDOFS = numCells
                     gradField = computeGradientFV2(varField, varConn, varCoord, varConStenDex)
                     
              # Store the gradients on target mesh
              try:
                     gradFileOut = nc_fid.createDimension(numDims, 3)
                     gradFileOut = nc_fid.createVariable(varGradientName, 'f8', (numDims, numDOFS))
                     gradFileOut[:] = gradField
              except Exception as exc:
                     print('Gradient variable already exists in ST field data file or: ', exc)
       
       nc_fid.close()
       endt = time.time()
       print('Time to compute/read gradients on target mesh (sec): ', endt - start)
       
       return gradField
       
if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison metrics!')
       print('Authors: Jorge Guerra, Paul Ullrich, 2019')

       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       varName, nc_fileSS, nc_fileS2T, nc_fileST, \
       mesh_fileS, mesh_fileT, \
       sourceMeshConfig, targetMeshConfig, sourceSE, targetSE = \
       parseCommandLine(sys.argv[1:])
       
       # Set the names for the auxiliary area and adjacency maps (NOT USER)
       varAreaName = 'cell_area'
       varJacoName = 'element_jacobians'
       varAdjaName = 'cell_edge_adjacency'
       varGradientName = 'FieldGradient'
              
       # Read in raw vertex/connectivity data from mesh files
       varCoordS, varConS, numEdgesS, numCellsS, numDimsS, numVertsS = \
       loadMeshData(mesh_fileS, sourceMeshConfig, sourceSE)

       varCoordT, varConT, numEdgesT, numCellsT, numDimsT, numVertsT = \
       loadMeshData(mesh_fileT, targetMeshConfig, targetSE)
              
       # Read in source and target cell areas
       areaS = loadMeshAreas(mesh_fileS, varAreaName)
       areaT = loadMeshAreas(mesh_fileT, varAreaName)
              
       # Read in source and target Jacobian weights
       jacobiansS = loadMeshJacobians(mesh_fileS, varJacoName, sourceSE)
       jacobiansT = loadMeshJacobians(mesh_fileT, varJacoName, targetSE)
       
       # Read in source and target adjacency maps
       varConStenDexS = loadMeshAdjacencyMap(mesh_fileS, varAdjaName)
       varConStenDexT = loadMeshAdjacencyMap(mesh_fileT, varAdjaName)
       
       # Read in field variable data
       varSS = loadField(nc_fileSS, varName)
       varST = loadField(nc_fileST, varName)
       varS2T = loadField(nc_fileS2T, varName)
       
       # Read in or compute the respective gradients on target mesh
       gradST = loadFieldGradient(nc_fileST, varGradientName, varST, varConT, varCoordT, varConStenDexT, jacobiansT, numCellsT, numDimsT, targetSE)
       gradS2T = loadFieldGradient(nc_fileS2T, varGradientName, varS2T, varConT, varCoordT, varConStenDexT, jacobiansT, numCellsT, numDimsT, targetSE)
       
       varsOnTM = [varST, varS2T]
       gradientsOnTM = [gradST, gradS2T]

       #%%
       start = time.time()
       print('Computing all metrics...')
       # Global conservation metric
       massS2T, massST, L_g = computeGlobalConservation(varConS, varConT, varSS, varS2T, varST, areaS, areaT, jacobiansS, jacobiansT, sourceSE, targetSE)
       # Locality measure (returns an array for each target DOF)
       #L_local = computeLocalityMetric(varS2T, varST, varConn, varCoord)
       # Standard Error norms (L_1, L_2, L_inf)
       L_1, L_2, L_inf = computeStandardNorms(varConT, varS2T, varST, areaT, jacobiansT, targetSE)
       # Global Extrema preservation
       Lmin, Lmax = computeGlobalExtremaMetrics(varS2T, varST)
       # Local Extrema preservation
       Lmin_inf, Lmax_inf = \
       computeLocalExtremaMetrics(varConStenDexT, varConT, varCoordT, varS2T, varST, targetSE)
       # Gradient preservation
       gradientsOnTM = [gradST, gradS2T]
       varsOnTM = [varST, varS2T]
       H1, H1_2 = computeGradientPreserveMetrics(varConT, gradientsOnTM, varsOnTM, areaT, jacobiansT, targetSE)
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
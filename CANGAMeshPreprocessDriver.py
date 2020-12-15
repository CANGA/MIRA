#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
NAME
    Mesh pre-processing for FieldGen and MetricsDriver
PURPOSE
    Verification steps for mesh data. If grid coordinates (Cartesian sphere center),
    and cell/element connectivity exist, then a verification flag is set to the
    mesh files. This is true for Exodus .g files. If the format is SCRIP where
    connectivity maps are not given, then compute them and set the verification
    flag. Also compute cell areas and/or adjacency maps and sets a flag to the
    mesh files accordingly.
PROGRAMMER(S)
    Jorge Guerra, Vijay Mehadevan, Paul Ullrich
REVISION HISTORY
    
REFERENCES
'''
#%%
import shutil
import time
import sys, os, getopt
import math as mt
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from computeCoordConFastSCRIP import computeCoordConFastSCRIP
from computeFastAdjacencyStencil import computeFastAdjacencyStencil
from computeCoordConnGLL import computeCoordConnGLL
from computeAreaIntegral import computeAreaIntegral, computeAreaIntegralWithGQ, getGaussNodesWeights
from computeAreaIntegralSE import computeAreaIntegralSE
import computeSphericalCartesianTransforms as sphcrt

import multiprocessing
from multiprocessing import Process
from itertools import repeat

NTASKS = 8

# Parse the command line
def parseCommandLine(argv):
       
       # Mesh information files
       sampleMesh = ''
       ExodusSingleConn = False
       ExodusMultiConn = False
       SCRIPwithoutConn = False
       SCRIPwithConn = False
       SpectralElement = False
       forceRecompute = False
       # Polynomial order for spectral elements
       seOrder = 4
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['mesh=', \
                                         'force', \
                                         'ExodusSingleConn', 'ExodusMultiConn', 'SCRIPwithoutConn', \
                                         'SCRIPwithConn', 'SpectralElement', 'seorder='])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAMeshPreprocessDriver.py', \
                    '--mesh <meshFile>', \
                    '--<meshConfiguration>', \
                    '--<makeSEGrid>', \
                    '--force', \
                    '--seorder <polyOrder>')
              sys.exit(2)
              
       for opt, arg in opts:
              # Request for usage help
              if opt == '-h':
                     print('Command line not properly set:', \
                           'CANGAMeshPreprocessDriver.py', \
                           '--mesh <meshFile>', \
                           '--<meshConfiguration>', \
                           '--<makeSEGrid>', \
                           '--seorder <polyOrder>')
                     sys.exit()
              elif opt == '--mesh':
                     sampleMesh = arg
              elif opt == '--force':
                     forceRecompute = True
              elif opt == '--ExodusSingleConn':
                     ExodusSingleConn = True
              elif opt == '--ExodusMultiConn':
                     ExodusMultiConn = True
              elif opt == '--SCRIPwithoutConn':
                     SCRIPwithoutConn = True
              elif opt == '--SCRIPwithConn':
                     SCRIPwithConn = True
              elif opt == '--SpectralElement':
                     SpectralElement = True
              elif opt == '--seorder':
                     if int(arg)%2 == 0 and int(arg) < 5:
                         seOrder = int(arg)
                     else:
                         sys.exit("[FATAL] Error in option passed for --seorder. SE order must be in [2, 4]")
                     
       # Check that only one configuration is chosen
       configs = [ExodusSingleConn, ExodusMultiConn, SCRIPwithoutConn, SCRIPwithConn]
       numConfigs = sum(bool(x) for x in configs)
       if numConfigs > 1:
              print('ONE mesh configuration option must be set!')
              print('None of the options are set.')
              sys.exit(2)
       
       return sampleMesh, ExodusSingleConn, ExodusMultiConn, SCRIPwithoutConn, \
              SCRIPwithConn, SpectralElement, seOrder, forceRecompute

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison mesh pre-processor!')
       print('Authors: Jorge Guerra, Vijay Mahadevan, Paul Ullrich, 2019')
       
       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       mesh_file, ExodusSingleConn, ExodusMultiConn, SCRIPwithoutConn, SCRIPwithConn, \
       SpectralElement, seOrder, forceRecompute = parseCommandLine(sys.argv[1:])

       # Set the names for the auxiliary area and adjacency maps (NOT USER)
       varAreaName = 'cell_area'
       varAdjaName = 'cell_edge_adjacency'
       
       # Set the name of the augmented mesh file and copy from original
       mfname, mfext = os.path.splitext(mesh_file)
       outFileName = "{0}_{1}{2}".format(mfname, "enhanced", mfext)
       shutil.copy(mesh_file, outFileName)
       
       #%% Mesh processing
       if ExodusSingleConn:
              numEdges = 'num_nod_per_el1'
              numCells = 'num_el_in_blk1'
              numDims = 'cart_dims'
              connCell = 'element_corners_id'
              coordCell = 'grid_corners_cart'
              numVerts = 'grid_corners_size'
              
              # Open the .g mesh files for reading
              m_fid = Dataset(mesh_file, 'a')
              
              start = time.time()
              # Get connectivity and coordinate arrays (check for multiple connectivity)
              varCon = m_fid.variables['connect1'][:]
              varCoord = m_fid.variables['coord'][:]
              
              try:   
                     print('Storing connectivity and coordinate arrays from Exodus mesh files.')
                     meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
                     meshFileOut = m_fid.createDimension(numDims, 3)
                     meshFileOut = m_fid.createVariable(connCell, 'i4', (numCells, numEdges))
                     meshFileOut[:] = varCon
                     meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts))
                     meshFileOut[:] = varCoord
                     
              except RuntimeError:
                     print('Cell connectivity and grid vertices exist in mesh data file.') 
              
              endt = time.time()
              print('Time to precompute EXODUS single connectivity mesh info (sec): ', endt - start)
              
       elif ExodusMultiConn:
              numElTypes = 'num_el_blk'
              numDims = 'cart_dims'
              connCell = 'element_corners_id'
              coordCell = 'grid_corners_cart'
              numVerts = 'grid_corners_size'
              
              # Open the .g mesh files for reading
              m_fid = Dataset(mesh_file, 'a')
              
              start = time.time()
              # Get connectivity and coordinate arrays
              varConnList = []
              numVertList = []
              numConnBlocks = len(m_fid.dimensions[numElTypes])
              for cc in range(numConnBlocks):
                     # Get this connectivity array (El X corners)
                     connName = 'connect' + str(cc+1)
                     thisConn = m_fid.variables[connName][:]
                     # Get the number of corners for this connectivity block
                     numVertList.append(thisConn.size[1]) # Column dimension of connectivity
                     # Append to the list of connectivity blocks
                     varConnList.append(m_fid.variables[connName][:])
                     
              # Get the maximum number of vertices
              maxVerts = np.amax(np.array(numVertList))
              # Loop over the blocks again and pad columns up to the max vertices
              for cc in range(numConnBlocks):
                     numVert2Pad = maxVerts - numVertList[cc]
                     
                     if numVert2Pad == 0:
                            continue
                     
                     # Pad with redundant last coord ID up to the max vertices
                     thisPadding = np.matlib.repmat(varConnList[cc][:,-1], 1, numVert2Pad)
                     varConnList[cc] = np.hstack((varConnList[cc], thisPadding))
                     
              # Vertical stack of the connectivity lists
              varConn = np.vstack(tuple(varConnList))
              print(varConn)
              varCoord = m_fid.variables['coord'][:]
              
              try:   
                     print('Storing connectivity and coordinate arrays from Exodus mesh files.')
                     numEdges = 'num_nod_per_el'
                     numCells = 'num_el_in_blk'
                     meshFileOut = m_fid.createDimension(numEdges, maxVerts)
                     meshFileOut = m_fid.createDimension(numCells, varConn.size[0])
                     meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
                     meshFileOut = m_fid.createDimension(numDims, 3)
                     meshFileOut = m_fid.createVariable(connCell, 'i4', (numCells, numEdges))
                     meshFileOut[:] = varCon
                     meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts))
                     meshFileOut[:] = varCoord
                     
              except RuntimeError:
                     print('Cell connectivity and grid vertices exist in mesh data file.') 
              
              endt = time.time()
              print('Time to precompute EXODUS multi-connectivity mesh info (sec): ', endt - start)
              
       elif SCRIPwithoutConn:
              numEdges = 'grid_corners'
              numCells = 'grid_size'
              numDims = 'cart_dims'
              connCell = 'element_corners_id'
              coordCell = 'grid_corners_cart'
              numVerts = 'grid_corners_size'
              
              # Open the .nc SCRIP files for reading
              m_fid = Dataset(mesh_file, 'a')
              
              # Get the list of available variables
              varList = m_fid.variables.keys()
              
              # Get RAW (no ID) connectivity and coordinate arrays
              conLon = m_fid.variables['grid_corner_lon'][:]
              conLat = m_fid.variables['grid_corner_lat'][:]
              
              # Convert to radians if necessary
              if m_fid.variables['grid_corner_lon'].units == 'degrees':
                     conLon *= mt.pi / 180.0
                     
              if m_fid.variables['grid_corner_lat'].units == 'degrees':
                     conLat *= mt.pi / 180.0
                   
              start = time.time()

              print('Computing connectivity and coordinate arrays from raw SCRIP without connectivity')
              # Make coordinate and connectivity from raw SCRIP data
              varCoordLL, varCon = computeCoordConFastSCRIP(conLon, conLat)
              
              # Convert coordinates from lat/lon to Cartesian
              varCoord = sphcrt.computeLL2Cart(varCoordLL[:,1:4])
              varCoord = varCoord.T
              
              try:   
                     print('Storing connectivity and coordinate arrays from raw SCRIP')
                     meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
                     meshFileOut = m_fid.createDimension(numDims, 3)
                     meshFileOut = m_fid.createVariable(connCell, 'i4', (numCells, numEdges))
                     meshFileOut[:] = varCon
                     meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts))
                     meshFileOut[:] = varCoord
                     
              except RuntimeError:
                     print('Cell connectivity and grid vertices exist in mesh data file.') 
              
              endt = time.time()
              print('Time to precompute SCRIP mesh info (sec): ', endt - start)
              
       elif SCRIPwithConn:
              numEdges = 'ncorners'
              numCells = 'ncells'
              numDims = 'cart_dims'
              connCell = 'element_corners'
              coordCell = 'grid_corners_cart'
              
              # Open the .nc SCRIP files for reading
              m_fid = Dataset(mesh_file)
              
              # Get the list of available variables
              varList = m_fid.variables.keys()
              
              # Get RAW coordinate arrays and connectivity
              conLon = m_fid.variables['lon'][:]
              conLat = m_fid.variables['lat'][:]
              varCon = m_fid.variables[connCell][:]
              varCon = varCon.T
              
              # Convert to radians if necessary
              if m_fid.variables['lon'].units == 'degrees_east':
                     conLon *= mt.pi / 180.0
                     
              if m_fid.variables['lat'].units == 'degrees_north':
                     conLat *= mt.pi / 180.0
                     
              start = time.time()
              
              print('Computing connectivity and coordinate arrays from raw SCRIP with connectivity')
              # Make coordinate from raw SCRIP data
              varCoordLL = np.zeros((len(conLon), 3))
              varCoordLL[:,0] = conLon
              varCoordLL[:,1] = conLat
              varCoordLL[:,2] = np.ones(len(conLon))
              
              # Convert coordinates from lat/lon to Cartesian
              varCoord = sphcrt.computeLL2Cart(varCoordLL)
              varCoord = varCoord.T
              
              try:
                     print('Storing connectivity and coordinate arrays from raw SCRIP')
                     meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
                     meshFileOut = m_fid.createDimension(numDims, 3)
                     meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts))
                     meshFileOut[:] = varCoord
                     
              except RuntimeError:
                     print('Cell connectivity and grid vertices exist in mesh data file.') 
              
              endt = time.time()
              print('Time to precompute SCRIP mesh info (sec): ', endt - start)
              
       # Close the mesh file
       m_fid.close()

       # Open the mesh file for new data
       m_fid = Dataset(outFileName, 'a')
              
       #%% Adjacency processing
       
       start = time.time()
       print('Computing adjacency maps...')
       try:
              tempVar = m_fid.variables[varAdjaName]
              varInFile = True
       except KeyError:
              varInFile = False

       if not varInFile or forceRecompute:
              try:
                     if not varInFile:
                            meshFileOut = m_fid.createVariable(varAdjaName, 'i4', (numCells, numEdges))
                     else:
                            meshFileOut = m_fid.variables[varAdjaName]

                     print('Adjacency data computed/written to mesh file for the first time...')
                     # Compute adjacency maps for both meshes (source stencil NOT needed)
                     edgeNodeMap, edgeNodeKDTree, varConStenDex = computeFastAdjacencyStencil(varCon)
                     # Get the starting index for the adjecency information in varConStenDexT
                     adex = np.size(varConStenDex,1) - np.size(varCon,1) 

                     meshFileOut[:] = varConStenDex[:,adex:]
              except RuntimeError:
                     print('Adjacency variable already exists in mesh data file')
       else:
              varConStenDex = m_fid.variables[varAdjaName]

       endt = time.time()
       print('Time to precompute adjacency maps (sec): ', endt - start)
       
       #%% Area processing for FV models
       start = time.time()
       print('Computing mesh areas...')
       try:
              tempVar = m_fid.variables[varAreaName]
              varInFile = True
       except KeyError:
              varInFile = False

       if not varInFile or forceRecompute:

              try:
                     NEL = len(varCon)
                     if not varInFile:
                            meshFileOut = m_fid.createVariable(varAreaName, 'f8', (numCells, ))
                            area = np.zeros((NEL))
                     else:
                            meshFileOut = m_fid.variables[varAreaName]
                            area = np.zeros(meshFileOut.shape)

                     GN, GW = getGaussNodesWeights(6)

                     # Loop over each cell and get cell average
                     pool = multiprocessing.Pool(processes=NTASKS)
                     results = pool.starmap(computeAreaIntegralWithGQ, zip(repeat(1.0), [varCoord[:, varCon[ii,:] - 1] for ii in range(NEL)], repeat(GN), repeat(GW), repeat(False), repeat(True)))
                     pool.close()
                     pool.join()
                     varAreas  = np.array(results, dtype='f8')[:, 1]

                     meshFileOut[:] = np.ravel(varAreas)
              except RuntimeError:
                     print('Source areas already exist in mesh data file.')
       else:
              varConStenDex = m_fid.variables[varAreaName]
              
       endt = time.time()
       print('Time to precompute cell areas (sec): ', endt - start)
       
       #%% Make global GLL connectivity and coordinate arrays from varCon and varCoord
       start = time.time()
       print('Computing new GLL mesh grids and connectivity...')
       
       # Check element topology (must be quads)
       NGC = varCon.shape[1]
       if ((NGC == 4) and (SpectralElement)):
              print('GLL global coordinates and connectivity computed/written to mesh file for the first time...')
              NEL = len(varCon)
       
              # Compute number of grids per element
              if seOrder == 2:
                     NGED = 2
                     NGEL = 4
              elif seOrder == 4:
                     NGED = 4
                     NGEL = 16
              else:
                     NGED = 4
                     NGEL = 16
                     print('Assuming 4th order Spectral Elements') 
              
              # Compute the new GLL global coordinates and connectivity (by edges)
              edgeNodeMapGLL, varCoordGLL, varConGLL = \
                     computeCoordConnGLL(NEL, NGED, NGEL, varCoord, varCon, edgeNodeMap, edgeNodeKDTree, seOrder)
                     
              try:
                     print('Storing GLL connectivity and coordinate arrays.')
                     numVertsGLL = 'grid_gll_size'
                     numNodesGLL = 'num_gll_per_el1'
                     connCellGLL = 'element_gll_conn'
                     coordCellGLL = 'grid_gll_cart'
                     
                     meshFileOut = m_fid.createDimension(numVertsGLL, np.size(varCoordGLL, 1))
                     meshFileOut = m_fid.createDimension(numNodesGLL, NGEL)
                     meshFileOut = m_fid.createVariable(connCellGLL, 'i4', (numCells, numNodesGLL))
                     meshFileOut[:] = varConGLL
                     meshFileOut = m_fid.createVariable(coordCellGLL, 'f8', (numDims, numVertsGLL))
                     meshFileOut[:] = varCoordGLL
                     
              except RuntimeError:
                     print('Cell connectivity and grid vertices exist in mesh data file.')
       else:
              print('GLL global grid will NOT be computed. Elements are NOT regular quadrilaterals or --SpectralElement option not set.')
       
       endt = time.time()
       print('Time to precompute GLL grid/connectivity (sec): ', endt - start)
       
       #%% Jacobian/area processing for SE models
       start = time.time()
       print('Computing mesh element Jacobian weights...')
       
       NGC = varCon.shape[1]
       if ((NGC == 4) and (SpectralElement)):
              print('Element Jacobians computed/written to mesh file for the first time...')
              # Precompute the area weights and then look them up in the integral below
              NEL = len(varConGLL)
              jacobians = np.zeros(varConGLL.shape)
              for ii in range(NEL):
                     cdex = varConGLL[ii,:] - 1
                     thisCell = varCoordGLL[:,cdex.astype(int)]
                     elArea, jacobians[ii,:] = computeAreaIntegralSE(thisCell, 4)       
              try:   
                            print('Storing GLL Jacobian arrays.')
                            jacobiansGLL = 'element_jacobians'
                            numNodesGLL = 'num_gll_per_el1'
                            
                            meshFileOut = m_fid.createVariable(jacobiansGLL, 'f8', (numCells, numNodesGLL))
                            meshFileOut[:] = jacobians
                            
              except RuntimeError:
                     print('Cell connectivity and grid vertices exist in mesh data file.')
       else:
              print('GLL global grid will NOT be computed. Elements are NOT regular quadrilaterals or --SpectralElement option not set.')
              
       endt = time.time()
       print('Time to precompute GLL element Jacobian weights (sec): ', endt - start)
       
       #%% Close out the file       
       m_fid.close()

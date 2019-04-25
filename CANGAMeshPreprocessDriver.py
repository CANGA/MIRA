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
import time
import sys, getopt
import math as mt
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from computeCoordConFastSCRIP import computeCoordConFastSCRIP
from computeAreaAverage import computeAreaAverage

def computeCart2LL(cellCoord):
       # Loop over each cell centroid, extract (lon, lat)
       NC = np.size(cellCoord, axis=0)
       varLonLat = np.zeros((NC, 2))
       for ii in range(NC):
              RO = np.linalg.norm(cellCoord[ii,:])
              psi = mt.asin(1.0 / RO * cellCoord[ii,2])
              lam = mt.atan2(-cellCoord[ii,0], -cellCoord[ii,1]) + mt.pi
              varLonLat[ii,:] = [lam, psi]
       
       # OUTPUT IS IN RADIANS       
       return varLonLat

def computeLL2Cart(cellCoord):
       # Loop over the Lon/Lat coordinate array, extract Cartesian coords
       # Input array is [lon, lat, radius]
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
       
       # Mesh information files
       sampleMesh = ''
       ExodusSingleConn = False
       SCRIPwithoutConn = False
       SCRIPwithConn = False
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['mesh=',
                                         'ExodusSingleConn', 'SCRIPwithoutConn', \
                                         'SCRIPwithConn'])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAMeshPreprocessDriver.py', \
                    '--mesh <meshFile>', \
                    '--<meshConfiguration>')
              sys.exit(2)
              
       for opt, arg in opts:
              # Request for usage help
              if opt == '-h':
                     print('Command line not properly set:', \
                           'CANGAMeshPreprocessDriver.py', \
                           '--mesh <meshFile>', \
                           '--<meshConfiguration>')
                     sys.exit()
              elif opt == '--mesh':
                     sampleMesh = arg
              elif opt == '--ExodusSingleConn':
                     ExodusSingleConn = True
              elif opt == '--SCRIPwithoutConn':
                     SCRIPwithoutConn = True
              elif opt == '--SCRIPwithConn':
                     SCRIPwithConn = True
                     
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
       
       return sampleMesh, ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison mesh pre-processor!')
       print('Authors: Jorge Guerra, Paul Ullrich, 2019')
       
       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       mesh_file, ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn = parseCommandLine(sys.argv[1:])
       
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
              
              # Get connectivity and coordinate arrays (check for multiple connectivity)
              varCon = m_fid.variables['connect1'][:]
              varCoord = m_fid.variables['coord'][:]
              
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
              try:
                     print('Reading connectivity and coordinate arrays from raw SCRIP')
                     varCon = m_fid.variables[connCell][:]
                     varCoord = m_fid.variables[coordCell][:]
              except KeyError:
                     print('Computing connectivity and coordinate arrays from raw SCRIP')
                     # Make coordinate and connectivity from raw SCRIP data
                     varCoordLL, varCon = computeCoordConFastSCRIP(conLon, conLat)
                     
                     # Convert coordinates from lat/lon to Cartesian
                     varCoord = computeLL2Cart(varCoordLL[:,1:4])
                     varCoord = varCoord.T
                     
                     try:   
                            print('Storing connectivity and coordinate arrays from raw SCRIP')
                            meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
                            meshFileOut = m_fid.createDimension(numDims, 3)
                            meshFileOut = m_fid.createVariable(connCell, 'i4', (numCells, numEdges, ))
                            meshFileOut[:] = varCon
                            meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts ))
                            meshFileOut[:] = varCoord
                            
                     except RuntimeError:
                            print('Cell connectivity and grid vertices exist in mesh data file.') 
              
              endt = time.time()
              print('Time to read/precompute SCRIP mesh info (sec): ', endt - start)
              
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
              try:
                     print('Reading connectivity and coordinate arrays from raw SCRIP')
                     varCon = m_fid.variables[connCell][:]
                     varCoord = m_fid.variables[coordCell][:]
              except KeyError:
                     print('Computing connectivity and coordinate arrays from raw SCRIP')
                     # Make coordinate from raw SCRIP data
                     varCoordLL = np.zeros((len(conLon), 3))
                     varCoordLL[:,0] = conLon
                     varCoordLL[:,1] = conLat
                     varCoordLL[:,2] = np.ones(len(conLon))
                     
                     # Convert coordinates from lat/lon to Cartesian
                     varCoord = computeLL2Cart(varCoordLL)
                     varCoord = varCoord.T
                     
                     try:   
                            print('Storing connectivity and coordinate arrays from raw SCRIP')
                            meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
                            meshFileOut = m_fid.createDimension(numDims, 3)
                            meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts ))
                            meshFileOut[:] = varCoord
                            
                     except RuntimeError:
                            print('Cell connectivity and grid vertices exist in mesh data file.') 
              
              endt = time.time()
              print('Time to read/precompute SCRIP mesh info (sec): ', endt - start)
              
       #%% Area and adjacency processing
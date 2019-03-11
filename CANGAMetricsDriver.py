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
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

# Bring in all the different metric modules
from computeGradient2 import computeGradient2
from computeCoordConSCRIP import computeCoordConSCRIP
from computeFastAdjacencyStencil import computeFastAdjacencyStencil
from computeGlobalConservation import computeGlobalConservation
#from computeLocalityMetric import computeLocalityMetric
from computeAreaWeight import computeAreaWeight
from computeStandardNorms import computeStandardNorms
from computeGlobalExtremaMetrics import computeGlobalExtremaMetrics
from computeLocalExtremaMetrics import computeLocalExtremaMetrics
from computeGradientPreserveMetrics import computeGradientPreserveMetrics

# Save grid variable (area, adjacency) back to mesh file
def saveNewMeshInfo(gridFile, var2store):
       # Put the two variables into netcdf file
       
       return

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
       AreaAdjacentyPrecomp = False
       
       # Mesh data configuration
       meshConfig = ''
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['ss=','s2t=','st=','sm=','tm=', \
                                         'ExodusSingleConn', 'SCRIPwithoutConn', 'AreaAdjacentyPrecomp'])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAMEtricsDriver.py', \
                    '-ss <sourceSampledFile>', \
                    '-s2t <remappedFile>', \
                    '-st <targetSampledFile>', \
                    '-sm <sourceMesh>', \
                    '-tm <targetMesh>', \
                    '--<meshConfiguration>', '--<doPrecompute>')
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
                           '--<meshConfiguration>', '--<doPrecompute>')
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
              elif opt == '--AreaAdjacentyPrecomp':
                     AreaAdjacentyPrecomp = True
                     print('WARNING: Will carry out EXPENSIVE Area and Adjacency map computations.')
                     
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
       elif (ExodusSingleConn == False) & (SCRIPwithoutConn == False):
              print('Expecting only ONE mesh configuration option!')
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
       
       print('Welcome to CANGA remapping intercomparison metrics!')              
       print('Mesh and Variable data must be in NETCDF format.')
       
       return varName, sourceSampledFile, targetSampledFile, \
              remappedFile, sourceMesh, targetMesh, \
              ExodusSingleConn, SCRIPwithoutConn

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison metrics!')
       print('When running in an IDE, comment out command line parsing: lines 146-147.')
       # Global parameters
       #kdleafs = 100
       #COINCIDENT_TOLERANCE = 1.0E-14

       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       #varName, nc_fileSS, nc_fileS2T, nc_fileST, mesh_fileS, mesh_fileT, \
       #ExodusSingleConn, SCRIPwithoutConn = parseCommandLine(sys.argv[1:])
       
       #""" SET INPUT HERE FOR DEVELOPMENT TESTING
       # Set the mesh configuration (mutually exclusive):
       # ExodusSingleConn -> DEFAULT BEST (DEGENERATE POLYGONS OK)
       # ExodusMultiConn -> NOT IMPLEMENTED (DEGENERATE POLYGONS OK)
       # SCRIPwithoutConn -> UNFORTUNATE SECOND
       # SCRIPwithConn -> NOT IMPLEMENTED (READING METIS MESH INFO PROBLEMATIC)
       ExodusSingleConn = True
       #ExodusMultiConn = False
       SCRIPwithoutConn = False
       #SCRIPwithConn = False
       
       # Set flag for precomputations
       AreaAdjacentyPrecomp = False
       
       # Set the name of the field variable in question (scalar)
       varName = 'Psi'
       
       # Set the names for the auxiliary area and adjacency maps
       varAreaName = 'cell_area'
       varAdjaName = 'cell_edge_adjacency'
       
       # Field sampled at the source (SS)
       nc_fileSS = 'testdata_CSne30_np4_3.nc'
       # Field mapped from source to target (S2T)
       nc_fileS2T = 'testdata_CSne30_2_RLL1deg_np4_3.nc'
       #nc_fileS2T = 'testdata_CSne30_2_ICO64_np4_3.nc'
       # Field sampled at the target (ST)
       nc_fileST = 'testdata_RLL1deg_np4_3.nc'
       #nc_fileST = 'testdata_ICO64_np4_3.nc'
       #""" 
       
       if ExodusSingleConn:
              numEdges = 'num_nod_per_el1'
              numCells = 'num_el_in_blk1'
              # Source Exodus .g file
              mesh_fileS = 'outCSne30.g'
              # Target Exodus .g file
              mesh_fileT = 'outRLL1deg.g'
              #mesh_fileT = 'outICO64.g'
              
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
              # Source SCRIP file
              mesh_fileS = 'Grids/ne30np4_pentagons.091226.nc'
              # Target SCRIP file
              mesh_fileT = 'Grids/ne30np4_latlon.091226.nc'
              
              # Open the .nc SCRIP files for reading
              m_fidS = Dataset(mesh_fileS, 'r')
              m_fidT = Dataset(mesh_fileT, 'r')
              
              # Get the list of available variables
              varListS = m_fidS.variables.keys()
              varListT = m_fidT.variables.keys()
              
              # Get RAW (no ID) connectivity and coordinate arrays
              lonS = m_fidS.variables['grid_corner_lon'][:]
              latS = m_fidS.variables['grid_corner_lat'][:]
              lonT = m_fidT.variables['grid_corner_lon'][:]
              latT = m_fidT.variables['grid_corner_lat'][:]
              
              # Make coordinate and connectivity from raw SCRIP data
              start = time.time()
              varCoordS, varConS = computeCoordConSCRIP(lonS, latS)
              varCoordT, varConT = computeCoordConSCRIP(lonT, latT)
              
              endt = time.time()
              print('Time to precompute SCRIP mesh info (sec): ', endt - start)
              
       m_fidS.close()
       m_fidT.close()
       
       if AreaAdjacentyPrecomp:
              #%%
              start = time.time()
              print('Computing adjacency maps...')
              # Compute adjacency maps for both meshes (source stencil NOT needed)
              edgeNodeMapT, edgeCellMapT, cleanEdgeCellMapT, varConStenDexT = computeFastAdjacencyStencil(varConT)
              # Store the adjacency map in the original grid netcdf file (target mesh)
              m_fidT = Dataset(mesh_fileT, 'a')
              meshFileOut = m_fidT.createVariable(varAdjaName, 'i4', (numCells, numEdges, ))
              meshFileOut[:] = varConStenDexT[:,4:]
              m_fidT.close()
                     
              endt = time.time()
              print('Time to precompute adjacency maps (sec): ', endt - start)
              #%%
              start = time.time()
              print('Computing source and target mesh areas...')
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
                     
              areaS = np.ravel(areaS)
              areaT = np.ravel(areaT)
              
              # Store the grid cell areas in the original netcdf file (source and target)
              m_fidS = Dataset(mesh_fileS, 'a')
              meshFileOut = m_fidS.createVariable(varAreaName, 'f8', (numCells, ))
              meshFileOut[:] = areaS
              m_fidS.close()
              m_fidT = Dataset(mesh_fileT, 'a')
              meshFileOut = m_fidT.createVariable(varAreaName, 'f8', (numCells, ))
              meshFileOut[:] = areaT
              m_fidT.close()
              
              endt = time.time()
              print('Time to precompute mesh areas (sec): ', endt - start)
       else:
              start = time.time()
              print('Reading cell areas and adjacency maps...')
              
              m_fidS = Dataset(mesh_fileS, 'r')
              m_fidT = Dataset(mesh_fileT, 'r')
              
              # Read in stored areas
              areaS = m_fidS.variables[varAreaName][:]
              areaT = m_fidT.variables[varAreaName][:]
              
              # Read in stored adjacency
              varConStenDexT = np.concatenate((varConT, m_fidT.variables[varAdjaName][:]), axis=1)
              
              endt = time.time()
              print('Time to read areas and adjacencies (sec): ', endt - start)
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
       print('Computing scalar gradients for target sampled and regridded fields...')
       # Precompute the gradient operator on regridded and sampled target data

       varsOnTM = [varST, varS2T]
       gradientsOnTM, cellCoordT = computeGradient2(varsOnTM, varCoordT, varConStenDexT, areaT)
       
       endt = time.time()
       print('Time to compute gradients on target mesh (sec): ', endt - start)
       
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
       print('Global conservation: ', L_g)
       print('Global L1 error:     ', L_1)
       print('Global L2 error:     ', L_2)
       print('Global Linf error:   ', L_inf)
       print('Global max error:    ', Lmax)
       print('Global min error:    ', Lmin)
       print('Local max L1 error:  ', Lmax_1)
       print('Local max L2 error:  ', Lmax_2)
       print('Local max Lm error:  ', Lmax_inf)
       print('Local min L1 error:  ', Lmin_1)
       print('Local min L2 error:  ', Lmin_2)
       print('Local min Lm error:  ', Lmin_inf)
       print('Gradient semi-norm:  ', H1_2)
       print('Gradient full-norm:  ', H1)       
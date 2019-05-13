#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
NAME
    Global field generator for remapping intercomparison
PURPOSE
    Reads 2 mesh data files (Exodus or SCRIP) and evaluates any one of, or
    combination of 3 fields (TPW, Cloud Fraction, Terrain) derived from
    Spherical Harmonic expansions of satellite global composite data.
PROGRAMMER(S)
    Jorge Guerra, Paul Ullrich
REVISION HISTORY
    
REFERENCES
'''    
#%%
import shutil
import time
import sys, getopt
import pyshtools
import math as mt
import numpy as np
import plotly as py
import plotly.figure_factory as FF
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from computeAreaIntegral import computeAreaIntegral

#%% Utility functions

def computeSpectrum(ND, lfPower, hfPower, degIntersect):
       psd = np.zeros(ND)
       # Compute power spectrum array from coefficients
       degs = np.arange(ND, dtype=float)
       #degs[0] = np.inf
       degs[0] = 1.0E-8
       for ii in range(ND):
              if degs[ii] < degIntersect:
                     psd[ii] = lfPower[0] * np.power(degs[ii], lfPower[1]) + lfPower[2]
              elif degs[ii] >= degIntersect:
                     psd[ii] = hfPower[0] * np.power(degs[ii], hfPower[1]) + hfPower[2]
       
       return degs, psd

def computeCentroids(varCon, varCoord):
       # Loop over cells and get centroid vectors
       NC = np.size(varCon, axis=0)
       NP = np.size(varCon, axis=1)
       cellCoord = np.zeros((NC, 3))
       for ii in range(NC):
              # Centroid by averaging corner sphere-center vectors
              centroid = np.mat([0.0, 0.0, 0.0])
              for pp in range(NP):
                     ndex = varCon[ii,pp] - 1
                     centroid += varCoord[:,ndex]
                     
              centroid *= 1.0 / NP
       
              # Renormalize the centroid vector
              RO = np.linalg.norm(centroid)
              centroid *= 1.0 / RO
              
              # Store this centroid
              cellCoord[ii,:] = centroid 
       
       return cellCoord

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

def computeCentroidsLL(conLon, conLat):
       # Loop over rows of the corner array and get centroid
       NC = np.size(conLon, axis=0)
       NP = np.size(conLon, axis=1)
       cellCoord = np.zeros((NC, 2))
       for ii in range(NC):
              # Centroid by averaging corner sphere-center vectors
              centroid = np.mat([0.0, 0.0])
              for pp in range(NP):
                     centroid += [conLon[ii,pp], conLat[ii,pp]]
                     
              centroid *= 1.0 / NP
              
              # Store this centroid
              cellCoord[ii,:] = centroid 
       
       return cellCoord

def computeCellAverage(clm, varCon, varCoord, order, avg):
       # Compute the number of cells and initialize
       NEL = np.size(varCon, 0)
       varSample = np.zeros(NEL,)
       
       # Loop over each cell and get cell average
       for ii in range(NEL):
              # Handle degeneracies with numpy.unique on connectivity
              cdex = np.unique(varCon[ii,:]) - 1
              thisCell = varCoord[:,cdex]
                            
              varSample[ii] = computeAreaIntegral(clm, thisCell, order, avg, False)
       
       return varSample

def computeRandomizedCoefficients(ND):
       # Initialize the coefficients array
       coeffs = np.zeros((2,ND,ND))
       
       # Set the random integer seed
       seed = 384
       
       # Loop over ND (number of degrees)
       for kk in range(ND):
              nrand = np.ones((2, kk+1))
              # Initialize random numbers with number of coefficients at this degree 
              if kk == 0:
                     rand = (1103515245 * seed + 25214903917 + 12345) % 2147483647
              # Loop over the coefficients at this degree
              for ll in range(0, kk+1):
                     nrand[0,ll] = rand
                     rand = (1103515245 * rand + 25214903917 + 12345) % 2147483647
                     nrand[1,ll] = rand
                     rand = (1103515245 * rand + 25214903917 + 12345) % 2147483647
         
              # Turn the random set into double
              nrand = np.multiply(nrand, 1.0 / 2147483647.0)
              
              # Set the coefficients at degree kk+1
              coeffs[:2,kk,:kk+1] = 2.0 * np.add(2.0 * nrand[:2,:], -1.0)
              
       return coeffs

# Parse the command line
def parseCommandLine(argv):
       
       # Mesh information files
       sampleMesh = ''
       ExodusSingleConn = False
       SCRIPwithoutConn = False
       SCRIPwithConn = False
       
       # Sampling order
       sampleCentroid = False
       sampleOrder = 4
       
       # SET WHICH FIELDS TO EVALUATE
       EvaluateAll = False
       EvaluateTPW = False # Total Precipitable Water
       EvaluateCFR = False # Global Cloud Fraction
       EvaluateTPO = False # Global topography
       
       # Number of modes used up to 768
       numModes = 128
       
       try:
              opts, args = getopt.getopt(argv, 'hv:', \
                                        ['pm=', 'so=', 'nm=', 'EvaluateAll', \
                                         'EvaluateTPW', 'EvaluateCFR', 'EvaluateTPO', \
                                         'ExodusSingleConn', 'SCRIPwithoutConn', \
                                         'SCRIPwithConn'])
       except getopt.GetoptError:
              print('Command line not properly set:', \
                    'CANGAFieldGenDriver.py', \
                    '--pm <sampleMeshFile>', \
                    '--so <sampleOrderInteger>', \
                    '--nm <numberSHModesMax768>', \
                    '--<evaluateAllFields>', \
                    '--<evaluateTotalPrecipWater>', \
                    '--<evaluateCloudFraction>', \
                    '--<evaluateGlobalTerrain>', \
                    '--<meshConfiguration>')
              sys.exit(2)
              
       for opt, arg in opts:
              # Request for usage help
              if opt == '-h':
                     print('Command line not properly set:', \
                           'CANGAFieldGenDriver.py', \
                           '--pm <sampleMeshFile>', \
                           '--so <sampleOrderInteger>', \
                           '--nm <numberSHModesMax768>', \
                           '--<evaluateAllFields>', \
                           '--<evaluateTotalPrecipWater>', \
                           '--<evaluateCloudFraction>', \
                           '--<evaluateGlobalTerrain>', \
                           '--<meshConfiguration>')
                     sys.exit()
              elif opt == '--pm':
                     sampleMesh = arg
              elif opt == '--so':
                     if int(arg) == 1:
                            sampleCentroid = True
                     else:
                            if int(arg)%2 == 0 and int(arg) < 7:
                                sampleOrder = int(arg)
                            else:
                                sys.exit("[FATAL] Error in option passed for --so. Sample order must be in [2, 4, 6]")
              elif opt == '--nm':
                     numModes = int(arg)
              elif opt == '--EvaluateAll':
                     EvaluateAll = True
              elif opt == '--EvaluateTPW':
                     EvaluateTPW = True
              elif opt == '--EvaluateCFR':
                     EvaluateCFR = True
              elif opt == '--EvaluateTPO':
                     EvaluateTPO = True
              elif opt == '--ExodusSingleConn':
                     ExodusSingleConn = True
              elif opt == '--SCRIPwithoutConn':
                     SCRIPwithoutConn = True
              elif opt == '--SCRIPwithConn':
                     SCRIPwithConn = True
                     
       # Check that the number of modes requested doesn't exceed 768
       if numModes > 768:
              print('Setting maximum number of expansion modes: 768.')
              numModes = 768
                     
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
       
       print('Welcome to CANGA remapping intercomparison metrics!')              
       print('Mesh and Variable data must be in NETCDF format.')
       
       return sampleMesh, numModes, \
              sampleCentroid, sampleOrder, \
              EvaluateAll, EvaluateTPW, EvaluateCFR, EvaluateTPO, \
              ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison field generator!')
       print('Authors: Jorge Guerra, Paul Ullrich, 2019')
       
       # Parse the commandline! COMMENT OUT TO RUN IN IDE
       mesh_file, ND, sampleCentroid, sampleOrder, \
       EvaluateAll, EvaluateTPW, EvaluateCFR, EvaluateTPO, \
       ExodusSingleConn, SCRIPwithoutConn, SCRIPwithConn = parseCommandLine(sys.argv[1:])
       
       # FV mode
       avg = True

       """ SET INPUT HERE FOR DEVELOPMENT TESTING
       ND = 48
       # Set the mesh configuration (mutually exclusive):
       # ExodusSingleConn -> DEFAULT BEST (DEGENERATE POLYGONS OK)
       # ExodusMultiConn -> NOT IMPLEMENTED (DEGENERATE POLYGONS OK)
       # SCRIPwithoutConn -> UNFORTUNATE SECOND
       # SCRIPwithConn -> NOT IMPLEMENTED (READING METIS MESH INFO PROBLEMATIC)
       ExodusSingleConn = False
       #ExodusMultiConn = False
       SCRIPwithoutConn = True
       SCRIPwithConn = False
       
       # Sampling Exodus .g file
       #mesh_file = 'outCSne30.g'
       #mesh_file = 'outRLL1deg.g'
       #mesh_file = 'outICO64.g'
       
       # Sampling SCRIP file
       mesh_file = 'Grids/ne30np4_pentagons.091226.nc'
       #mesh_file = 'Grids/ne30np4_latlon.091226.nc'
       
       sampleCentroid = False
       sampleOrder = 6
       
       # SET WHICH FIELDS TO EVALUATE
       EvaluateAll = False
       EvaluateTPW = True # Total Precipitable Water
       EvaluateCFR = False # Global Cloud Fraction
       EvaluateTPO = False # Global topography
       """
       
       # Set the name for the new data file
       stripDir = mesh_file.split('/')
       onlyFilename = stripDir[len(stripDir)-1]
       data_file = 'testdata_' + (onlyFilename.split('.'))[0]
       print('New data will be stored in (prefix): ', data_file)
       
       print('Number of SH degrees for sampling set to: ', ND)
       
       if ExodusSingleConn:
              # Open the .g mesh files for reading
0              m_fid = Dataset(mesh_file)
              
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
                   
              start = time.time()
              try:
                     print('Reading connectivity and coordinate arrays from raw SCRIP')
                     varCon = m_fid.variables[connCell][:]
                     varCoord = m_fid.variables[coordCell][:]
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
              m_fid = Dataset(mesh_file, 'a')
              
              # Get the list of available variables
              varList = m_fid.variables.keys()
              
              # Get RAW (no ID) connectivity and coordinate arrays
              varCon = m_fid.variables[connCell][:]
              varCon = varCon.T
              
              start = time.time()
              try:
                     print('Reading coordinate arrays from raw SCRIP')
                     varCoord = m_fid.variables[coordCell][:]
              except:
                     print('PRE-PROCESSING NOT DONE ON THIS MESH FILE!')
              
              endt = time.time()
              print('Time to read SCRIP mesh info (sec): ', endt - start)
                            
       # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
       # Compute Centroids
       varCent = computeCentroids(varCon, varCoord)
       
       # Compute Lon/Lat coordinates from centroids
       varLonLat = computeCart2LL(varCent)
       
       # Convert to degrees from radians
       varLonLat_deg = 180.0 / mt.pi * varLonLat
       varLonLat_deg = 180.0 / mt.pi * varLonLat
       
       m_fid.close()
              
       #%% Begin the SH reconstructions
       if EvaluateTPW or EvaluateAll:
              start = time.time()
              print('Computing Total Precipitable Water on sampling mesh...')
              # Set the power spectrum coefficients
              lfPower = [5.84729561e+04, -2.91678103e-04, -5.83966265e+04]
              hfPower = [2.17936330e+02, -1.99788552e+00, -7.94469251e-04]
              degIntersect = 1.8161917668847762
              # Compute the parent power spectrum for TPW
              degsTPW, psdTPW = computeSpectrum(ND, lfPower, hfPower, degIntersect)

              # Set the low degree coefficients (large scale structures)
              coeffsLD_TPW = np.array([[2.45709150e+01, 0.0, 0.0, 0.0], \
                              [4.00222122e+00, 2.39412571e+00, 0.0, 0.0], \
                              [-1.36433589e+01, 3.90520866e-03, 4.70350344e-01, 0.0], \
                              [-3.54931720e+00, -1.23629157e+00, 4.01454924e-01, 1.76782768e+00]])
       
              # Initialize SHCoeffs with a randomized realization of coefficients
              clmTPW = pyshtools.SHCoeffs.from_random(psdTPW, seed=384)

              # Compute the randomized coefficients and update instance of SHCoeffs
              clmTPW.coeffs = computeRandomizedCoefficients(ND)
              
              # Force the coefficients to have the same power as the given spectrum
              power_per_l = pyshtools.spectralanalysis.spectrum(clmTPW.coeffs, normalization='4pi', unit='per_l')
              clmTPW.coeffs *= np.sqrt(psdTPW[0:ND] / power_per_l)[np.newaxis, :, np.newaxis]
              
              # Combine the coefficients, low degree from data and high degree randomized
              clmTPW.coeffs[0,0:4,0:4] = coeffsLD_TPW
              
              # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
              # Expand the coefficients and check the field
              if sampleCentroid:              
                     TPWvar = clmTPW.expand(lon=varLonLat_deg[:,0], lat=varLonLat_deg[:,1])
              else:
                     TPWvar = computeCellAverage(clmTPW, varCon, varCoord, sampleOrder, avg)
              
              # Compute rescaled data from 0.0 to max
              minTPW = np.amin(TPWvar)
              maxTPW = np.amax(TPWvar)
              deltaTPW = abs(maxTPW - minTPW)
              TPWvar = np.add(TPWvar, -minTPW)
              TPWvar *= maxTPW / deltaTPW
              endt = time.time()
              print('Time to compute TPW (mm): ', endt - start)
       #%%                     
       if EvaluateCFR or EvaluateAll:
              start = time.time()
              print('Computing Cloud Fraction on sampling mesh...')
              # Set the power spectrum coefficients
              lfPower = [8.38954430e+00, -1.85962382e-04, -8.38439294e+00]
              hfPower = [1.25594628e-01, -1.99203168e+00,  1.91763519e-06]
              degIntersect = 8.322269484619733
              # Compute the parent power spectrum for CFR
              degsCFR, psdCFR = computeSpectrum(ND, lfPower, hfPower, degIntersect)
              
              # Set the low degree coefficients (large scale structures)
              coeffsLD_CFR = np.array([[6.65795054e-01, 0.0, 0.0, 0.0], \
                              [-2.45480409e-02, 2.24697424e-02, 0.0, 0.0], \
                              [5.72322008e-02, 3.41184683e-02, -7.71082815e-03, 0.0], \
                              [1.86562455e-02, 4.34697733e-04, 8.91735978e-03, -5.53756958e-03]])
       
              # Initialize SHCoeffs with a randomized realization of coefficients
              clmCFR = pyshtools.SHCoeffs.from_random(psdCFR, seed=384)
              
              # Compute the randomized coefficients and update instance of SHCoeffs
              clmCFR.coeffs = computeRandomizedCoefficients(ND)
              
              # Force the coefficients to have the same power as the given spectrum
              power_per_l = pyshtools.spectralanalysis.spectrum(clmCFR.coeffs, normalization='4pi', unit='per_l')
              clmCFR.coeffs *= np.sqrt(psdCFR[0:ND] / power_per_l)[np.newaxis, :, np.newaxis]
              
              # Combine the coefficients, low degree from data and high degree randomized
              clmCFR.coeffs[0,0:4,0:4] = coeffsLD_CFR
              
              # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
              # Expand the coefficients and check the field
              if sampleCentroid:              
                     CFRvar = clmCFR.expand(lon=varLonLat_deg[:,0], lat=varLonLat_deg[:,1])
              else:
                     CFRvar = computeCellAverage(clmCFR, varCon, varCoord, sampleOrder, avg)
 
              # Compute rescaled data from 0.0 to max
              minCFR = np.amin(CFRvar)
              maxCFR = np.amax(CFRvar)
              deltaCFR = abs(maxCFR - minCFR)
              CFRvar = np.add(CFRvar, -minCFR)
              CFRvar *= maxCFR / deltaCFR
              #  Set all values greater than 1.0 to 1.0 (creates discontinuities)
              CFRvar[CFRvar >= 1.0] = 1.0
              
              endt = time.time()
              print('Time to compute CFR (0.0 to 1.0): ', endt - start)
       #%%       
       if EvaluateTPO or EvaluateAll:
              start = time.time()
              print('Computing Terrain on sampling mesh...')
              # Set the power spectrum coefficients
              lfPower = [1.79242815e+05, -4.28193211e+01,  7.68040558e+05]
              hfPower = [9.56198160e+06, -1.85485966e+00, -2.63553217e+01]
              degIntersect = 3.8942282772035255
              # Compute the parent power spectrum for CFR
              degsTPO, psdTPO = computeSpectrum(ND, lfPower, hfPower, degIntersect)
              
              # Set the low degree coefficients (large scale structures)
              coeffsLD_TPO = np.array([[-2.38452711e+03, 0.0, 0.0, 0.0], \
                              [-6.47223253e+02, -6.06453097e+02, 0.0, 0.0], \
                              [5.67394318e+02, 3.32672611e+02, -4.17639577e+02, 0.0], \
                              [1.57403492e+02, 1.52896988e+02, 4.47106726e+02, -1.40553447e+02]])
                         
              # Initialize SHCoeffs with a randomized realization of coefficients
              clmTPO = pyshtools.SHCoeffs.from_random(psdTPO, seed=384)
              
              # Compute the randomized coefficients and update instance of SHCoeffs
              clmTPO.coeffs = computeRandomizedCoefficients(ND)
              
              # Force the coefficients to have the same power as the given spectrum
              power_per_l = pyshtools.spectralanalysis.spectrum(clmTPO.coeffs, normalization='4pi', unit='per_l')
              clmTPO.coeffs *= np.sqrt(psdTPO[0:ND] / power_per_l)[np.newaxis, :, np.newaxis]
              
              # Combine the coefficients, low degree from data and high degree randomized
              clmTPO.coeffs[0,0:4,0:4] = coeffsLD_TPO
              
              # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
              # Expand the coefficients and check the field
              if sampleCentroid:              
                     TPOvar = clmTPO.expand(lon=varLonLat_deg[:,0], lat=varLonLat_deg[:,1])
              else:
                     TPOvar = computeCellAverage(clmTPO, varCon, varCoord, sampleOrder, avg)
              
              # Rescale to -1.0 to 1.0
              minTPO = np.amin(TPOvar)
              maxTPO = np.amax(TPOvar)
              deltaTPO = abs(maxTPO - minTPO)
              TPOvar = np.add(TPOvar, -0.5 * (maxTPO + minTPO))
              TPOvar *= 2.0 / deltaTPO
              
              # Rescale topography to real Earth max/min
              minTPO = -10994.0 # Depth at Challenger Deep
              maxTPO = 8848.0 # Elevation of Mt. Everest ASL
              deltaTPO = abs(maxTPO - minTPO)
              TPOvar *= (0.5 * deltaTPO)
              TPOvar += 0.5 * (maxTPO + minTPO)
              
              endt = time.time()
              print('Time to compute TPO (m): ', endt - start)
              
       #%% Copy grid files and store the new test data (source and target)
       outFileName = data_file
       if EvaluateAll:
              outFileName = outFileName + '_TPW_CFR_TPO.nc'
       elif EvaluateTPW:
              outFileName = outFileName + '_TPW.nc'
       elif EvaluateCFR:
              outFileName = outFileName + '_CFR.nc'
       elif EvaluateTPO:
              outFileName = outFileName + '_TPO.nc'
       else:
              outFileName = outFileName + '.nc'
              
       shutil.copy(mesh_file, outFileName)
       
       # write lon, lat, and test data variables
       data_fid = Dataset(outFileName, 'a')
       
       # Set the dimension name depending on the mesh file format
       if ExodusSingleConn:
              numCells = 'num_el_in_blk1'
       elif SCRIPwithoutConn:
              numCells = 'grid_size'
       elif SCRIPwithConn:
              numCells = 'ncells'
              
       # Process the sampling file
       if SCRIPwithConn:
              lonNC = data_fid.createVariable('nlon', 'f8', (numCells,))
              lonNC[:] = varLonLat_deg[:,0]
              latNC = data_fid.createVariable('nlat', 'f8', (numCells,))
              latNC[:] = varLonLat_deg[:,1]
       else:
              lonNC = data_fid.createVariable('lon', 'f8', (numCells,))
              lonNC[:] = varLonLat_deg[:,0]
              latNC = data_fid.createVariable('lat', 'f8', (numCells,))
              latNC[:] = varLonLat_deg[:,1]
       
       if EvaluateTPW or EvaluateAll:
              TPWNC = data_fid.createVariable('TotalPrecipWater', 'f8', (numCells,))
              TPWNC[:] = TPWvar
       if EvaluateCFR or EvaluateAll:
              CFRNC = data_fid.createVariable('CloudFraction', 'f8', (numCells,))
              CFRNC[:] = CFRvar
       if EvaluateTPO or EvaluateAll:
              TPONC = data_fid.createVariable('Topography', 'f8', (numCells,))
              TPONC[:] = TPOvar
       
       # Close the files out.
       data_fid.close()

       '''

       #%% Check the data with triangular surface plot
       points2D = varLonLat
       tri = Delaunay(points2D)
       simplices = tri.simplices       
       
       #%% Plot Total Precipitable Water
       if EvaluateTPW or EvaluateAll:
              fig1 = FF.create_trisurf(x=varLonLat[:,0], y=varLonLat[:,1], z=TPWvar, height=800, width=1200, \
                                       simplices=simplices, colormap="Portland", plot_edges=False, \
                                       title="Total Precipitable Water Check (mm)", aspectratio=dict(x=1, y=1, z=0.3))
              py.offline.plot(fig1, filename='TPW' + data_file + '.html')
       #%% Plot Cloud Fraction
       if EvaluateCFR or EvaluateAll:
              fig1 = FF.create_trisurf(x=varLonLat[:,0], y=varLonLat[:,1], z=CFRvar, height=800, width=1200, \
                                       simplices=simplices, colormap="Portland", plot_edges=False, \
                                       title="Cloud Fraction Check (0.0-1.0)", aspectratio=dict(x=1, y=1, z=0.3))
              py.offline.plot(fig1, filename='CFR' + data_file + '.html')
       #%% Plot Topography
       if EvaluateTPO or EvaluateAll:
              fig1 = FF.create_trisurf(x=varLonLat[:,0], y=varLonLat[:,1], z=TPOvar, height=800, width=1200, \
                                       simplices=simplices, colormap="Portland", plot_edges=False, \
                                       title="Global Topography (m)", aspectratio=dict(x=1, y=1, z=0.3))
              py.offline.plot(fig1, filename='TPO' + data_file + '.html')
       '''
       #%% Check the evaluated spectra
       '''
       fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(12, 10), tight_layout=True)
       # Plot the TPW spectrum
       newPSD = pyshtools.spectralanalysis.spectrum(clmTPW.coeffs, unit='per_l')
       ax0.plot(degsTPW, psdTPW, 'k')
       ax0.plot(degsTPW, newPSD, 'r--')
       ax0.set_title('Total Precipitable Water - Evaluated PSD')
       ax0.set(yscale='log', xscale='log', ylabel='Power')
       ax0.grid(b=True, which='both', axis='both')
       # Plot the Cloud Fraction spectrum
       newPSD = pyshtools.spectralanalysis.spectrum(clmCFR.coeffs, unit='per_l')
       ax1.plot(degsCFR, psdCFR, 'k')
       ax1.plot(degsCFR, newPSD, 'r--')
       ax1.set_title('Global Cloud Fraction - Evaluated PSD')
       ax1.set(yscale='log', xscale='log', ylabel='Power')
       ax1.grid(b=True, which='both', axis='both')
       # Plot the Topography spectrum
       newPSD = pyshtools.spectralanalysis.spectrum(clmTPO.coeffs, unit='per_l')
       ax2.plot(degsTPO, psdTPO, 'k')
       ax2.plot(degsTPO, newPSD, 'r--')
       ax2.set_title('Global Topography Data - Evaluated PSD')
       ax2.set(yscale='log', xscale='log', xlabel='Spherical harmonic degree', ylabel='Power')
       ax2.grid(b=True, which='both', axis='both')
       plt.show()
       '''      

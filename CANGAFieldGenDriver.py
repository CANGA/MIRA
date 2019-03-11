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
import pyshtools
import math as mt
import numpy as np
import plotly as py
import plotly.figure_factory as FF
import plotly.graph_objs as go
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

#%% Utility functions

def computeSpectrum(ND, lfPower, hfPower, degIntersect):
       psd = np.zeros(ND)
       # Compute power spectrum array from coefficients
       degs = np.arange(ND, dtype=float)
       #degs[0] = np.inf
       degs[0] = 1.0E-6
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
              
       return varLonLat

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

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison field generator!')
       print('When running in an IDE, comment out command line parsing: lines 146-147.')
       
       ND = 256
       print('Number of SH degrees for sampling set to: ', ND)

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
       
       # SET WHICH FIELDS TO EVALUATE ON BOTH MESHES
       EvaluateAll = True
       EvaluateTPW = False # Total Precipitable Water
       EvaluateCFR = False # Global Cloud Fraction
       EvaluateTPO = False # Global topography
       
       if ExodusSingleConn:
              # Source Exodus .g file
              mesh_fileS = 'outCSne30.g'
              # Target Exodus .g file
              mesh_fileT = 'outRLL1deg.g'
              #mesh_fileT = 'outICO64.g'
              
              # Set a file name for new test data
              data_fileS = 'testdata_' + (mesh_fileS.split('.'))[0]
              data_fileT = 'testdata_' + (mesh_fileT.split('.'))[0]
              
              # Open the .g mesh files for reading
              g_fidS = Dataset(mesh_fileS)
              g_fidT = Dataset(mesh_fileT)
              
              # Get connectivity and coordinate arrays (check for multiple connectivity)
              varConS = g_fidS.variables['connect1'][:]
              varCoordS = g_fidS.variables['coord'][:]
              varConT = g_fidT.variables['connect1'][:]
              varCoordT = g_fidT.variables['coord'][:]
              
              # Compute Centroids
              varCentS = computeCentroids(varConS, varCoordS)
              varCentT = computeCentroids(varConT, varCoordT)
              
              # Compute Lon/Lat coordinates from centroids
              varLonLatS = computeCart2LL(varCentS)
              varLonLatT = computeCart2LL(varCentT)
              
              # Convert to degrees from radians
              varLonLatS_deg = 180.0 / mt.pi * varLonLatS
              varLonLatT_deg = 180.0 / mt.pi * varLonLatT
              
              g_fidS.close()
              g_fidT.close()
              
       elif SCRIPwithoutConn:
              # Source SCRIP file
              mesh_fileS = 'Grids/ne30np4_pentagons.091226.nc'
              # Target SCRIP file
              mesh_fileT = 'Grids/ne30np4_latlon.091226.nc'
              
              # Set a file name for new test data
              data_fileS = 'testdata_' + (mesh_fileS.split('.'))[0]
              data_fileT = 'testdata_' + (mesh_fileT.split('.'))[0]
              
              # Open the .nc SCRIP files for reading
              s_fidS = Dataset(mesh_fileS)
              s_fidT = Dataset(mesh_fileT)
              
              # Get the list of available variables
              varListS = s_fidS.variables.keys()
              varListT = s_fidT.variables.keys()
              
              # Get RAW (no ID) connectivity and coordinate arrays
              conLonS = s_fidS.variables['grid_corner_lon'][:]
              conLatS = s_fidS.variables['grid_corner_lat'][:]
              conLonT = s_fidT.variables['grid_corner_lon'][:]
              conLatT = s_fidT.variables['grid_corner_lat'][:]
              
              # Compute centroids from Lat/Lon corners
              varLonLatS = computeCentroidsLL(conLonS, conLatS)
              varLonLatT = computeCentroidsLL(conLonT, conLatT)
              
              # Convert to degrees from radians
              varLonLatS_deg = 180.0 / mt.pi * varLonLatS
              varLonLatT_deg = 180.0 / mt.pi * varLonLatT
              
              s_fidS.close()
              s_fidT.close()
              
       #%% Begin the SH reconstructions
       if EvaluateTPW or EvaluateAll:
              start = time.time()
              print('Computing Total Precipitable Water on source and target...')
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
       
              # Compute a randomized realization of coefficients
              clmTPW = pyshtools.SHCoeffs.from_random(psdTPW, exact_power=True, seed=512)
              
              # Combine the coefficients, low degree from data and high degree randomized
              clmTPW.coeffs[0,0:4,0:4] = coeffsLD_TPW
              
              # Expand the coefficients and check the field              
              TPWvarS = clmTPW.expand(lon=varLonLatS_deg[:,0], lat=varLonLatS_deg[:,1])
              TPWvarT = clmTPW.expand(lon=varLonLatT_deg[:,0], lat=varLonLatT_deg[:,1])
              # Compute rescaled data from 0.0 to max
              minTPW = np.amin(TPWvarS)
              maxTPW = np.amax(TPWvarS)
              deltaTPW = abs(maxTPW - minTPW)
              TPWvarS = np.add(TPWvarS, -minTPW)
              TPWvarS *= maxTPW / deltaTPW
              minTPW = np.amin(TPWvarT)
              maxTPW = np.amax(TPWvarT)
              deltaTPW = abs(maxTPW - minTPW)
              TPWvarT = np.add(TPWvarT, -minTPW)
              TPWvarT *= maxTPW / deltaTPW
              endt = time.time()
              print('Time to compute TPW (mm): ', endt - start)
                            
       if EvaluateCFR or EvaluateAll:
              start = time.time()
              print('Computing Cloud Fraction on source and target...')
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
       
              # Compute a randomized realization of coefficients
              clmCFR = pyshtools.SHCoeffs.from_random(psdCFR, exact_power=True, seed=512)
              
              # Combine the coefficients, low degree from data and high degree randomized
              clmCFR.coeffs[0,0:4,0:4] = coeffsLD_CFR
              
              # Expand the coefficients and check the field              
              CFRvarS = clmCFR.expand(lon=varLonLatS_deg[:,0], lat=varLonLatS_deg[:,1])
              CFRvarT = clmCFR.expand(lon=varLonLatT_deg[:,0], lat=varLonLatT_deg[:,1])
              # Compute rescaled data from 0.0 to max
              minCFR = np.amin(CFRvarS)
              maxCFR = np.amax(CFRvarS)
              deltaCFR = abs(maxCFR - minCFR)
              CFRvarS = np.add(CFRvarS, -minCFR)
              CFRvarS *= maxCFR / deltaCFR
              minCFR = np.amin(CFRvarT)
              maxCFR = np.amax(CFRvarT)
              deltaCFR = abs(maxCFR - minCFR)
              CFRvarT = np.add(CFRvarT, -minCFR)
              CFRvarT *= maxCFR / deltaCFR
              #  Set all values greater than 1.0 to 1.0 (creates discontinuities)
              CFRvarS[CFRvarS >= 1.0] = 1.0
              CFRvarT[CFRvarT >= 1.0] = 1.0
              
              endt = time.time()
              print('Time to compute CFR (0.0 to 1.0): ', endt - start)
              
       if EvaluateTPO or EvaluateAll:
              start = time.time()
              print('Computing Terrain on source and target...')
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
                         
              # Compute a randomized realization of coefficients
              clmTPO = pyshtools.SHCoeffs.from_random(psdTPO, exact_power=True, seed=512)
              
              # Combine the coefficients, low degree from data and high degree randomized
              clmTPO.coeffs[0,0:4,0:4] = coeffsLD_TPO
              
              # Expand the coefficients and check the field              
              TPOvarS = clmTPO.expand(lon=varLonLatS_deg[:,0], lat=varLonLatS_deg[:,1])
              TPOvarT = clmTPO.expand(lon=varLonLatT_deg[:,0], lat=varLonLatT_deg[:,1])
              # DO NOT rescale the topography
              endt = time.time()
              print('Time to compute TPO (m): ', endt - start)
              
       #%% Copy grid files and store the new test data (source and target)
       outFileNameS = data_fileS
       outFileNameT = data_fileT
       if EvaluateAll:
              outFileNameS = outFileNameS + '_TPW_CFR_TPO.nc'
              outFileNameT = outFileNameT + '_TPW_CFR_TPO.nc'
       elif EvaluateTPW:
              outFileNameS = outFileNameS + '_TPW.nc'
              outFileNameT = outFileNameT + '_TPW.nc'
       elif EvaluateCFR:
              outFileNameS = outFileNameS + '_CFR.nc'
              outFileNameT = outFileNameT + '_CFR.nc'
       elif EvaluateTPO:
              outFileNameS = outFileNameS + '_TPO.nc'
              outFileNameT = outFileNameT + '_TPO.nc'
       else:
              outFileNameS = outFileNameS + '.nc'
              outFileNameT = outFileNameT + '.nc'
              
       shutil.copy(mesh_fileS, outFileNameS)
       shutil.copy(mesh_fileT, outFileNameT)
       
       # write lon, lat, and test data variables
       data_fidS = Dataset(outFileNameS, 'a')
       data_fidT = Dataset(outFileNameT, 'a')
       
       # Set the dimension name depending on the mesh file format
       if ExodusSingleConn:
              numCells = 'num_el_in_blk1'
       elif SCRIPwithoutConn:
              numCells = 'grid_size'
              
       # Process the source file
       lonNC = data_fidS.createVariable('lon', 'f8', (numCells,))
       lonNC[:] = varLonLatS_deg[:,0]
       latNC = data_fidS.createVariable('lat', 'f8', (numCells,))
       latNC[:] = varLonLatS_deg[:,1]
       # Process the target file
       lonNC = data_fidT.createVariable('lon', 'f8', (numCells,))
       lonNC[:] = varLonLatT_deg[:,0]
       latNC = data_fidT.createVariable('lat', 'f8', (numCells,))
       latNC[:] = varLonLatT_deg[:,1]
       
       if EvaluateTPW or EvaluateAll:
              TPWNC = data_fidS.createVariable('TotalPrecipWater', 'f8', (numCells,))
              TPWNC[:] = TPWvarS
              TPWNC = data_fidT.createVariable('TotalPrecipWater', 'f8', (numCells,))
              TPWNC[:] = TPWvarT
       if EvaluateCFR or EvaluateAll:
              CFRNC = data_fidS.createVariable('CloudFraction', 'f8', (numCells,))
              CFRNC[:] = CFRvarS
              CFRNC = data_fidT.createVariable('CloudFraction', 'f8', (numCells,))
              CFRNC[:] = CFRvarT
       if EvaluateTPO or EvaluateAll:
              TPONC = data_fidS.createVariable('Topography', 'f8', (numCells,))
              TPONC[:] = TPOvarS
              TPONC = data_fidT.createVariable('Topography', 'f8', (numCells,))
              TPONC[:] = TPOvarT
       
       # Close the files out.
       data_fidS.close()
       data_fidT.close()

       #%% Check the data with triangular surface plot
       points2D = varLonLatT
       tri = Delaunay(points2D)
       simplices = tri.simplices       
       
       fig1 = FF.create_trisurf(x=varLonLatT[:,0], y=varLonLatT[:,1], z=TPWvarT, height=800, width=1200, \
                                simplices=simplices, colormap="Portland", plot_edges=False, \
                                title="Total Precipitable Water Check (mm)", aspectratio=dict(x=1, y=1, z=0.3))
       py.offline.plot(fig1, filename='TPW' + (mesh_fileT.split('.'))[0] + '.html')
       fig1 = FF.create_trisurf(x=varLonLatT[:,0], y=varLonLatT[:,1], z=CFRvarT, height=800, width=1200, \
                                simplices=simplices, colormap="Portland", plot_edges=False, \
                                title="Cloud Fraction Check (0.0-1.0)", aspectratio=dict(x=1, y=1, z=0.3))
       py.offline.plot(fig1, filename='CFR' + (mesh_fileT.split('.'))[0] + '.html')
       fig1 = FF.create_trisurf(x=varLonLatT[:,0], y=varLonLatT[:,1], z=TPOvarT, height=800, width=1200, \
                                simplices=simplices, colormap="Portland", plot_edges=False, \
                                title="Global Topography (m)", aspectratio=dict(x=1, y=1, z=0.3))
       py.offline.plot(fig1, filename='TPO' + (mesh_fileT.split('.'))[0] + '.html')
       
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
              
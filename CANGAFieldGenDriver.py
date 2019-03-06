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
import sys, getopt
import time
import pyshtools
import math as mt
import numpy as np
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
              lam = mt.atan(cellCoord[ii,0] / cellCoord[ii,1])
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
       
       ND = 501
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
              
       elif SCRIPwithoutConn:
              # Source SCRIP file
              mesh_fileS = 'Grids/ne30np4_pentagons.091226.nc'
              # Target SCRIP file
              mesh_fileT = 'Grids/ne30np4_latlon.091226.nc'
              
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
              
       #%% Begin the SH reconstructions
       if EvaluateTPW or EvaluateAll:
              # Set the power spectrum coefficients
              lfPower = [5.84729561e+04, -2.91678103e-04, -5.83966265e+04]
              hfPower = [2.17936330e+02, -1.99788552e+00, -7.94469251e-04]
              degIntersect = 1.8161917668847762
              # Compute the parent power spectrum for TPW
              degsTPW, psdTPW = computeSpectrum(ND, lfPower, hfPower, degIntersect)
              # "Fix" the end value
              psdTPW[0] = 2.0 * psdTPW[1]              
              # Compute a randomized realization of coefficients
              clmTPW = pyshtools.SHCoeffs.from_random(psdTPW, exact_power=True, seed=512)
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
                            
       if EvaluateCFR or EvaluateAll:
              # Set the power spectrum coefficients
              lfPower = [8.38954430e+00, -1.85962382e-04, -8.38439294e+00]
              hfPower = [1.25594628e-01, -1.99203168e+00,  1.91763519e-06]
              degIntersect = 8.322269484619733
              # Compute the parent power spectrum for CFR
              degsCFR, psdCFR = computeSpectrum(ND, lfPower, hfPower, degIntersect)
              # "Fix" the end value
              psdCFR[0] = 2.0 * psdCFR[1]
              # Compute a randomized realization of coefficients
              clmCFR = pyshtools.SHCoeffs.from_random(psdCFR, exact_power=True, seed=512)
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
              
       if EvaluateTPO or EvaluateAll:
              # Set the power spectrum coefficients
              lfPower = [1.79242815e+05, -4.28193211e+01,  7.68040558e+05]
              hfPower = [9.56198160e+06, -1.85485966e+00, -2.63553217e+01]
              degIntersect = 3.8942282772035255
              # Compute the parent power spectrum for CFR
              degsTPO, psdTPO = computeSpectrum(ND, lfPower, hfPower, degIntersect)
              # "Fix" the end value
              psdTPO[0] = 2.0 * psdTPO[1]              
              # Compute a randomized realization of coefficients
              clmTPO = pyshtools.SHCoeffs.from_random(psdTPO, exact_power=True, seed=512)
              # Expand the coefficients and check the field              
              TPOvarS = clmTPO.expand(lon=varLonLatS_deg[:,0], lat=varLonLatS_deg[:,1])
              TPOvarT = clmTPO.expand(lon=varLonLatT_deg[:,0], lat=varLonLatT_deg[:,1])
              # Compute rescaled data from 0.0 to max
              minTPO = np.amin(TPOvarS)
              maxTPO = np.amax(TPOvarS)
              deltaTPO = abs(maxTPO - minTPO)
              TPOvarS = np.add(TPOvarS.data, -minTPO)
              TPOvarS *= maxTPO / deltaTPO
              minTPO = np.amin(TPOvarT)
              maxTPO = np.amax(TPOvarT)
              deltaTPO = abs(maxTPO - minTPO)
              TPOvarT = np.add(TPOvarT, -minTPO)
              TPOvarT *= maxTPO / deltaTPO
              
       #%% Copy grid files and store the new test data (source and target)
       #fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(12, 10), tight_layout=True)
              
       #%% Check the evaluated spectra
       '''
       fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, figsize=(12, 10), tight_layout=True)
       # Plot the TPW spectrum
       ax0.plot(degsTPW, psdTPW, 'k')
       ax0.set_title('Total Precipitable Water - Evaluated PSD')
       ax0.set(yscale='log', xscale='log', ylabel='Power')
       ax0.grid(b=True, which='both', axis='both')
       # Plot the Cloud Fraction spectrum
       ax1.plot(degsCFR, psdCFR, 'k')
       ax1.set_title('Global Cloud Fraction - Evaluated PSD')
       ax1.set(yscale='log', xscale='log', ylabel='Power')
       ax1.grid(b=True, which='both', axis='both')
       # Plot the Topography spectrum
       ax2.plot(degsTPO, psdTPO, 'k')
       ax2.set_title('Global Topography Data - Evaluated PSD')
       ax2.set(yscale='log', xscale='log', xlabel='Spherical harmonic degree', ylabel='Power')
       ax2.grid(b=True, which='both', axis='both')
       plt.show()
       '''
              
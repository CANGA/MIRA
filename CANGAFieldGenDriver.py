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
import numpy as np
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/

#%% Utility functions

def computeSpectrumTPW():
       PSD = 0.0
       
       return PSD

def computeSpectrumCFR():
       PSD = 0.0
       
       return PSD

def computeSpectrumTPO():
       PSD = 0.0
       
       return PSD

def computeCentroids(varCon, varCoord):
       # Loop over cells and get centroid vectors
       
       return

def computeCentroidsRLL(lon, lat):
       # Loop over rows of the corner array and get centroid
       
       return

if __name__ == '__main__':
       print('Welcome to CANGA remapping intercomparison field generator!')
       print('When running in an IDE, comment out command line parsing: lines 146-147.')

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
              lonS, latS = computeCart2RLL(varCentS)
              lonT, latT = computeCart2RLL(varCentT)
              
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
              lonS = s_fidS.variables['grid_corner_lon'][:]
              latS = s_fidS.variables['grid_corner_lat'][:]
              lonT = s_fidT.variables['grid_corner_lon'][:]
              latT = s_fidT.variables['grid_corner_lat'][:]
              
              # Compute centroids from Lat/Lon corners
              lonS, latS = computeCentroidsRLL(lonS, latS)
              lonT, latT = computeCentroidsRLL(lonT, latT)
              
       #%% Begin the SH reconstructions
       if EvaluateTPW or EvaluateAll:
              # Set the power spectrum coefficients
              
       if EvaluateCFR or EvaluateAll:
              # Set the power spectrum coefficients
              
       if EvaluateTPO or EvaluateAll:
              # Set the power spectrum coefficients
              
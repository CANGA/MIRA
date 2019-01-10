#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 15:50:24 2019

 Test of the grid and connectivity from SCRIP data routine. Uses MPAS grids.

@author: TempestGuerra
"""
from netCDF4 import Dataset
from computeCoordConSCRIP import computeCoordConSCRIP

# SCRIP netcdf file
scp_file = 'Grids/ne30np4_pentagons.091226.nc'

# Open the .nc SCRIP files for reading
s_fid = Dataset(scp_file)

# Get the list of available variables
varList = s_fid.variables.keys()

# Get RAW (no ID) connectivity and coordinate arrays
lon = s_fid.variables['grid_corner_lon'][:]
lat = s_fid.variables['grid_corner_lat'][:]

# Make coordinate and connectivity from raw SCRIP data
# 2D lat/lon grid on a circle
varCoord, varCon = computeCoordConSCRIP(lon, lat)

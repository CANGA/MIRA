# CANGA-Metrics-Driver
Python driver for the CANGA-ROO project. Developed with Python 3.6 using Spyder IDE

Main file: CANGAMetricsDriver.py

Initial development test data is provided with this repo in the form of 3 .nc files and 2 .g (mesh files).

REQUIRES: http://code.google.com/p/netcdf4-python/ Python NetCDF IO modules

MAIN ASSUMPTIONS:
1) Non-coincident grid array from mesh file
2) Non-coincident element connectivity array from mesh file (single type 3 or 4 node)
    ** Support for multiple connectivity arrays coming soon **
3) Main variable is defined at cell centers only
4) Gradient metric is implemented using a linear reconstruction between adjacent cells
5) Locality metric is NOT yet implemented since it requires the regridding operator explicitly

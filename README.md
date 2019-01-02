# CANGA-Metrics-Driver
Python driver for the CANGA-ROO project. Developed with Python 3.6 using Spyder IDE

Main file: CANGAMetricsDriver.py

TO TEST: User change the following lines in the main function of the file above

       # Set the name of the field variable in question (scalar)
       varName = 'Psi'
       
       # Field sampled at the source (SS)
       nc_fileSS = 'testdata_CSne30_np4_3.nc'
       # Field mapped from source to target (S2T)
       nc_fileS2T = 'testdata_CSne30_2_RLL1deg_np4_3.nc'
       # Field sampled at the target (ST)
       nc_fileST = 'testdata_RLL1deg_np4_3.nc'
       
       # Source Exodus .g file (only needed for global conservation)
       exo_fileS = 'outCSne30.g'
       # Target Exodus .g file
       exo_fileT = 'outRLL1deg.g'

Initial development test data is provided with this repo in the form of 3 .nc files and 2 .g (mesh files).

REQUIRES: http://code.google.com/p/netcdf4-python/ Python NetCDF IO modules

MAIN ASSUMPTIONS:
1) Non-coincident grid array from mesh file
2) Non-coincident element connectivity array from mesh file (single type 3 or 4 node)
    ** Support for multiple connectivity arrays coming soon **
3) Main variable is defined at cell centers only
4) Gradient metric is implemented using a linear reconstruction between adjacent cells
5) Locality metric is NOT yet implemented since it requires the regridding operator explicitly

TO DO:
1) Grid array and connectivity generation from raw coordinate lists in the SCRIP format
2) Support for fields on GLL grids

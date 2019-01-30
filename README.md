# CANGA-Metrics-Driver
Python driver for the CANGA-ROO project. Developed with Python 3.6 using Spyder IDE

Main file: CANGAMetricsDriver.py

TO TEST: Run the following command line with the included data

1) Finite volume remap from cube-sphere to RLL grid, analytical test 3, using Exodus files

# python CANGAMetricsDriver.py -v Psi --ss testdata_CSne30_np4_3.nc --s2t testdata_CSne30_2_RLL1deg_np4_3.nc --st testdata_RLL1deg_np4_3.nc --sm outCSne30.g --tm outRLL1deg.g --ExodusSingleConn

2) Finite volume remap from cube-sphere to icosahedral grid, analytical test 3, using Exodus files

# python CANGAMetricsDriver.py -v Psi --ss testdata_CSne30_np4_3.nc --s2t testdata_CSne30_2_ICO64_np4_3.nc --st testdata_ICO64_np4_3.nc --sm outCSne30.g --tm outICO64.g --ExodusSingleConn

NOTES:
- variable name as given in data netcdf file
- last argument giving the configuration is necessary

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
1) Test support for SCRIP grid data
2) Support for fields on GLL grids

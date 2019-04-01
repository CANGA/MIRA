# CANGA-Metrics-Driver
Python driver for the CANGA-ROO project. Developed with Python 3.6 using Spyder IDE

Main field generator file: CANGAFieldGenerator.py
Main metrics file: CANGAMetricsDriver.py

TESTING SEQUENCE *EXAMPLE*:
1) Mesh files:
   - Source mesh file eg. "outCSne30.g"
   - Target mesh file et. "outRLL1deg.g"
   - BOTH MESH FILES MUST BE MATCHING FORMAT: EXODUS OR SCRIP

2) To generate sampled data on source and target:
   - python CANGAFieldGenDriver.py --pm outCSne30.g --so 6 --nm 768 --EvaluateAll --ExodusSingleConn
   - python CANGAFieldGenDriver.py --pm outRLL1deg.g --so 6 --nm 768 --EvaluateAll --ExodusSingleConn
   - In this example the sampling order is maximum 6th order and number of modes is maximum 768 which is roughly equivalent      to 0.25 degree resolution. This is the maximum supported based on expansions of satellite data that yield ~1000+ modes      for Total Precipitable Water, Cloud Fraction, and Topography.
   - This will generate NEW files with the prefix 'testdata_' and suffix 'TPW_CFR_TPO' if --EvaluateAll is set. Individual        fields can also be set with --EvaluateTPW, --EvaluateCFR, --EvaluateTPO.
   - The following NETCDF files are created: testdata_outCSne30_TPW_CFR_TPO.nc AND testdata_outRLL1deg_TPW_CFR_TPO.nc
   - The new data files are augmented copies of the original mesh data files keeping all metadata consistent.
   
3) Remap testdata files generated in 2) using TempestRemap
   - Run the following shell scripts provided (modify as you see fit): CSne30_2_RLL1deg_Remap_TPW_CFR_TPO.sh
   - STOP: This step assumes that you have downloaded and compiled tempest remap from:                                            https://github.com/ClimateGlobalChange/tempestremap
   - YOUR FAVORITE REMAPPER WILL LIKELY DO SOMETHING DIFFERENT ENOUGH TO BREAK THINGS WITH REGARD TO NETCDF DETAILS SO            FURTHER TESTING/FIXING WILL BE REQUIRED TO ACCOUNT FOR MORE DIVERSE USER CASES.
   - "Anywho..." as Paul U. might say, the script above will generate the following: CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc
  
4) To generate metrics on remapped data:
   - You now have the 3 essential files: 
   ** FIELD SAMPLED ON THE SOURCE MESH: testdata_outCSne30_TPW_CFR_TPO.nc
   ** FIELD SAMPLED ON THE TARGET MESH: testdata_outRLL1deg_TPW_CFR_TPO.nc
   ** FIELD REMAPPED SOURCE TO TARGET: CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc
   - Run the following shell script to generate metrics: CSne30_2_RLL1deg_np4_TPW_CFR_TPO.nc
   - IF this is the first run, the code will go into a LONG precomputation of cell areas, adjacencies, and field gradients        to be stored in the mesh file (areas and adjacency) and data file (gradients).
   - IF the code finds precomputed data, execuation will be MUCH faster. 

5) Output for TPW (AFTER precomputations) should look like this (on a 2015 MacBook Pro...):

Welcome to CANGA remapping intercomparison metrics!
Authors: Jorge Guerra, Paul Ullrich, 2019
Computing/reading adjacency maps...
Time to precompute adjacency maps (sec):  0.002249002456665039
Computing source and target mesh areas...
Time to precompute/read mesh areas (sec):  0.0022902488708496094
Time to read NC and Exodus data (sec):  0.0025649070739746094
Computing or reading gradients for target sampled and regridded fields...
Time to compute/read gradients on target mesh (sec):  0.0046977996826171875
Computing all metrics...
Time to execute metrics (sec):  41.048941135406494
Global conservation: -1.342620427744820e-11
Global L1 error:     -5.375076262147747e-02
Global L2 error:     5.818910116017171e-02
Global Linf error:   9.500633281658483e-02
Global max error:    -2.221532280056399e-02
Global min error:    8.323725405955121e-03
Local max L1 error:  1.422438389080596e-02
Local max L2 error:  3.435138652959414e-02
Local max Lm error:  1.506976424475193e-01
Local min L1 error:  -9.498159999055940e-05
Local min L2 error:  1.239117251460247e-03
Local min Lm error:  2.771707334760692e-02
Gradient semi-norm:  1.891627764536130e+00
Gradient full-norm:  1.892522541756834e+00

NOTES:
- variable name as given in data netcdf file, check with "ncdump -c <filename>".
- last argument giving the configuration is necessary.
- SCRIP mesh files provided by Phil J. have a variety of dimension/variable definitions, all of which, need to be included
  for this work.

REQUIRES: 
1) http://code.google.com/p/netcdf4-python/ Python NetCDF IO modules
2) https://shtools.oca.eu/shtools/ Python spherical harmonic tools package
3) Numpy
4) Scipy (KDTree search)
5) https://plot.ly/python/ Plotly (Fancy, web-based plotting)

MAIN ASSUMPTIONS:
1) Non-coincident grid array from mesh file
2) Non-coincident, convex, element connectivity array from mesh file (single type 3 or 4 node)
    ** Support for multiple connectivity arrays coming soon **
3) Main variable is defined at cell centers only (FV to FV mappings)
4) Gradient metric is implemented using a linear reconstruction between adjacent cells
5) Locality metric requires 2 remapping applications to measure on the same (source) grid

TO DO:
1) Test support for SCRIP grid data (provided by Phil J.)
2) Support for fields on GLL grids

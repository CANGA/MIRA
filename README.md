# CANGA-Metrics-Driver
Python driver for the CANGA-ROO project. Developed with Python 3.6 using Spyder IDE

- Main mesh preprocessing file: CNAGAMeshPreprocessDriver.py
- Main field generator file: CANGAFieldGenerator.py
- Main metrics file: CANGAMetricsDriver.py

TESTING SEQUENCE *EXAMPLE*:
1) Mesh files:
   - Source mesh file eg. "outCSne60.g"
   - Target mesh file et. "outRLL1deg.g"
   - THE TWO MESH FILES NEED NOT BE FORMAT CONSISTENT. THEY WILL BE PREPROCESSED INDEPENDENTLY.
   
2) Preprocess source and target meshes:
   - python CANGAMeshPreprocessDriver.py --mesh outCSne60.g --ExodusSingleConn --SpectralElement
   - python CANGAMeshPreprocessDriver.py --mesh outRLL1deg.g --ExodusSingleConn
   - This will write the two variables varCoord (Cartesian coordinates of global nodes where the field is evaluated i.e. cell centers in the FV case or GLL nodes in the FE case) and varCon (Global connectivity list for cells/elements) to the respective mesh file.
   - Output will be:
   - Augmented source mesh file: "outCSne60.g_RIPPED"
   - Augmented target mesh file: "outRLL1deg.g_RIPPED"
   - The --SpectralElement option only applies to quadrilaterals (code will check and abort otherwise) and generates augmented global DOF coordinates and connectivity. It also generates the global Jacobian weights for each global DOF. This facilitates integrals through the susequent metrics computations.

3) To generate sampled data on source and target:
   - python CANGAFieldGenDriver.py --pm outCSne60.g_RIPPED --so 6 --nm 512 --EvaluateAll --ExodusSingleConn
   - THIS WILL OUTPUT ON CELL CENTERS USING A 6TH ORDER TRIANGULAR QUADRATURE TO COMPUTE ELEMENT AVERAGE.
   - python CANGAFieldGenDriver.py --pm outCSne60.g_RIPPED --so 6 --nm 512 --EvaluateAll --ExodusSingleConn --SpectralElement
   - THIS WILL OUTPUT ON GLOBAL GLL NODES WITH NO QUADRATURE APPLIED SO --so OPTION IS IGNORED.
   - python CANGAFieldGenDriver.py --pm outRLL1deg.g_RIPPED --so 6 --nm 512 --EvaluateAll --ExodusSingleConn
   - In this example the sampling order is maximum 6th order and number of modes is maximum 768 which is roughly equivalent      to 0.25 degree resolution. This is the maximum supported based on expansions of satellite data that yield ~1000+ modes      for Total Precipitable Water, Cloud Fraction, and Topography.
   - This will generate NEW files with the prefix 'testdata_' and suffix 'TPW_CFR_TPO' if --EvaluateAll is set. Individual        fields can also be set with --EvaluateTPW, --EvaluateCFR, --EvaluateTPO.
   - The following NETCDF files are created: testdata_outCSne60_TPW_CFR_TPO.nc AND testdata_outRLL1deg_TPW_CFR_TPO.nc
   - The new data files are augmented copies of the original mesh data files keeping all metadata consistent.
   
4) Remap testdata files generated in 2) using TempestRemap
   - Run the following shell scripts provided (modify as you see fit): CSne60_2_RLL1deg_Remap_TPW_CFR_TPO.sh
   - STOP: This step assumes that you have downloaded and compiled tempest remap from:                                            https://github.com/ClimateGlobalChange/tempestremap
   - YOUR FAVORITE REMAPPER WILL LIKELY DO SOMETHING DIFFERENT ENOUGH TO BREAK THINGS WITH REGARD TO NETCDF DETAILS SO            FURTHER TESTING/FIXING WILL BE REQUIRED TO ACCOUNT FOR MORE DIVERSE USER CASES.
   - "Anywho..." as Paul U. might say, the script above will generate the following: CSne60_2_RLL1deg_np4_TPW_CFR_TPO.nc
  
5) To generate metrics on remapped data:
   - You now have the 3 essential files: 
   ** FIELD SAMPLED ON THE SOURCE MESH: testdata_outCSne60_TPW_CFR_TPO.nc
   ** FIELD SAMPLED ON THE TARGET MESH: testdata_outRLL1deg_TPW_CFR_TPO.nc
   ** FIELD REMAPPED SOURCE TO TARGET: CSne60_2_RLL1deg_np4_TPW_CFR_TPO.nc
   - Run the following shell script to generate metrics: CSne60_2_RLL1deg_np4_TPW_CFR_TPO.nc
   - IF this is the first run, the code will go into a LONG precomputation of cell areas, adjacencies, and field gradients        to be stored in the mesh file (areas and adjacency) and data file (gradients).
   - IF the code finds precomputed data, execuation will be MUCH faster. 

6) Output for TPW (AFTER precomputations) should look like this (on a 2015 MacBook Pro...):


NOTES:
- variable name as given in data netcdf file, check with "ncdump -c <filename>".
- last argument giving the configuration is necessary.
- SCRIP mesh files provided by Phil J. have a variety of dimension/variable definitions, all of which, need to be included
  for this work.

REQUIRES: 
1) http://code.google.com/p/netcdf4-python/ Python NetCDF IO modules: "pip install netcdf4"
2) https://shtools.oca.eu/shtools/ Python spherical harmonic tools package: "pip install pyshtools"
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

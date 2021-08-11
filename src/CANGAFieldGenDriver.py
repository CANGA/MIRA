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
    Jorge Guerra, Vijay Mehadevan, Paul Ullrich
REVISION HISTORY

REFERENCES
'''
# %%
import shutil
import time
import sys
import getopt
import pyshtools
import math as mt
import numpy as np
from numpy import matlib
import plotly as py
import plotly.figure_factory as FF
from scipy.spatial import Delaunay
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
from computeAreaIntegral import computeAreaIntegral, computeAreaIntegralWithGQ, getGaussNodesWeights
import computeSphericalCartesianTransforms as sphcrt

import multiprocessing
from multiprocessing import Process
from itertools import repeat


# %% Utility functions

def computeSpectrum(ND, lfPower, hfPower, degIntersect):
    psd = np.zeros(ND)
    # Compute power spectrum array from coefficients (Power Law assumed)
    degs = np.arange(ND, dtype=float)
    # degs[0] = np.inf
    degs[0] = 1.0E-8

    # Check that we aren't fitting a constant function (Terrain)
    for ii in range(ND):
        if degs[ii] < degIntersect:
            if lfPower[1] > -5.0:
                psd[ii] = lfPower[0] * np.power(degs[ii], lfPower[1]) + lfPower[2]
            else:
                psd[ii] = lfPower[2]
        elif degs[ii] >= degIntersect:
            if hfPower[1] > -5.0:
                psd[ii] = hfPower[0] * np.power(degs[ii], hfPower[1]) + hfPower[2]
            else:
                psd[ii] = hfPower[2]

    return degs, psd


def evaluate_field_a2(lon, lat):
    # thisVar = (2.0 + np.cos(dFLonLat[1]) * np.cos(dFLonLat[1]) * np.cos(2.0 * dFLonLat[0])) # test == 1
    # thisVar = (2.0 + (np.sin(2.0 * dFLonLat[1]))**16.0 * np.cos(16.0 * dFLonLat[0])) # test == 2
    # print(lon, lat, (2.0 + np.cos(lat) * np.cos(lat) * np.cos(2.0 * lon)))
    return (2.0 + np.cos(lat) * np.cos(lat) * np.cos(2.0 * lon))


def computeCellAverageSerial(clm, varCon, varCoord, order, avg):
    # Compute the number of cells and initialize
    NEL = np.size(varCon, 0)
    varSample = np.zeros(NEL)

    # Loop over each cell and get cell average
    for ii in range(NEL):
        # NP.UNIQUE SORTS AND DESTROYS CONNECTIVITY CELL NORMALS!!!
        cdex = varCon[ii, :] - 1
        thisCell = varCoord[:, cdex]

        varSample[ii] = computeAreaIntegral(clm, thisCell, order, avg, False)

    return varSample


def computeCellAverage(clm, varCon, varCoord, order, avg, nprocs):

    # return computeCellAverageSerial(clm, varCon, varCoord, order, avg)

    # Compute the number of cells and initialize
    NEL = np.size(varCon, 0)
    varSample = np.zeros(NEL,)

    GN, GW = getGaussNodesWeights(order)

    # Loop over each cell and get cell average
    pool = multiprocessing.Pool(processes=nprocs)
    results = pool.starmap(computeAreaIntegralWithGQ, zip(
        repeat(clm), [varCoord[:, varCon[ii, :] - 1] for ii in range(NEL)], repeat(GN), repeat(GW), repeat(avg), repeat(False)))
    pool.close()
    pool.join()
    varSample = np.array(results, dtype='f8')[:, 0]
    varAreas = np.array(results, dtype='f8')[:, 1]

    return varSample


def computeRandomizedCoefficients(ND):
    # Initialize the coefficients array
    coeffs = np.zeros((2, ND, ND))

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
            nrand[0, ll] = rand
            rand = (1103515245 * rand + 25214903917 + 12345) % 2147483647
            nrand[1, ll] = rand
            rand = (1103515245 * rand + 25214903917 + 12345) % 2147483647

        # Turn the random set into double
        nrand = np.multiply(nrand, 1.0 / 2147483647.0)

        # Set the coefficients at degree kk+1
        coeffs[:2, kk, :kk+1] = 2.0 * np.add(2.0 * nrand[:2, :], -1.0)

    return coeffs


def computeNormalizedCoefficients(N, psd, coeffsLD):
    # Initialize SHCoeffs with a randomized realization of coefficients
    clm = pyshtools.SHCoeffs.from_random(psd, seed=384)

    # Compute the randomized coefficients and update instance of SHCoeffs
    clm.coeffs = computeRandomizedCoefficients(ND)

    # Force the coefficients to have the same power as the given spectrum
    power_per_l = pyshtools.spectralanalysis.spectrum(clm.coeffs, normalization='4pi', unit='per_l')
    clm.coeffs *= np.sqrt(psd[0:ND] * np.reciprocal(power_per_l))[np.newaxis, :, np.newaxis]

    # Combine the coefficients, low degree from data and high degree randomized
    clm.coeffs[0, 0:4, 0:4] = coeffsLD

    # Returns the SH coefficients object
    return clm

# Parse the command line


def parseCommandLine(argv):

    # Mesh information files
    sampleMesh = ''
    ExodusSingleConn = False
    ExodusMultiConn = False
    SCRIPwithoutConn = False
    SCRIPwithConn = False
    SpectralElement = False

    # Sampling order
    sampleCentroid = False
    sampleOrder = 4

    # SET WHICH FIELDS TO EVALUATE
    EvaluateAll = False
    EvaluateTPW = False  # Total Precipitable Water
    EvaluateCFR = False  # Global Cloud Fraction
    EvaluateTPO = False  # Global topography
    EvaluateA1 = False  # Analytical function 1
    EvaluateA2 = False  # Analytical function 2

    ShowPlots = False  # Whether we want to show the profile plots for variables

    # Number of modes used up to 512
    numModes = 32

    # Pseudo-random number generator seed
    seed = 384

    # Number of processes to use for sampling
    nprocs = 1

    def usage():
        print('Driver Usage:\n',
              'CANGAFieldGenDriver.py',
              '--pm <sampleMeshFile>',
              '--so <sampleOrderInteger>',
              '--nm <numberSHModesMax768>',
              '--rseed <randnumSeed>',
              '--evaluateAllFields',
              '--evaluateTotalPrecipWater',
              '--evaluateCloudFraction',
              '--evaluateGlobalTerrain',
              '--evaluateA1',
              '--evaluateA2',
              '--showPlots',
              '--meshConfiguration',
              '--SpectralElementMesh',
              '--processes <nprocs>')

    try:
        opts, args = getopt.getopt(argv, 'hv:',
                                   ['pm=', 'so=', 'nm=', 'rseed=', 'evaluateAllFields',
                                    'evaluateTotalPrecipWater', 'evaluateCloudFraction', 'evaluateGlobalTerrain',
                                    'evaluateA1', 'evaluateA2', 'showPlots',
                                    'ExodusSingleConn', 'ExodusMultiConn', 'SCRIPwithoutConn',
                                    'SCRIPwithConn', 'SpectralElementMesh', 'processes='])
    except getopt.GetoptError:
        print('Command line arguments were not properly set or error in parsing.\n')
        usage()
        sys.exit(2)

    for opt, arg in opts:
        # Request for usage help
        if opt == '-h':
            usage()
            sys.exit()
        elif opt == '--pm':
            sampleMesh = arg
        elif opt == '--so':
            if int(arg) == 1:
                sampleOrder = int(arg)
                sampleCentroid = True
            else:
                if int(arg) % 2 == 0 and int(arg) < 200:
                    sampleOrder = int(arg)
                else:
                    sys.exit("[FATAL] Error in option passed for --so. Sample order must be \in (0, 200)")
        elif opt == '--nm':
            numModes = int(arg)
        elif opt == '--rseed':
            seed = int(arg)
        elif opt == '--evaluateAllFields':
            EvaluateAll = True
        elif opt == '--evaluateTotalPrecipWater':
            EvaluateTPW = True
        elif opt == '--evaluateCloudFraction':
            EvaluateCFR = True
        elif opt == '--evaluateGlobalTerrain':
            EvaluateTPO = True
        elif opt == '--evaluateA1':
            EvaluateA1 = True
        elif opt == '--evaluateA2':
            EvaluateA2 = True
        elif opt == '--ExodusSingleConn':
            ExodusSingleConn = True
        elif opt == '--ExodusMultiConn':
            ExodusMultiConn = True
        elif opt == '--SCRIPwithoutConn':
            SCRIPwithoutConn = True
        elif opt == '--SCRIPwithConn':
            SCRIPwithConn = True
        elif opt == '--SpectralElementMesh':
            SpectralElement = True
        elif opt == '--showPlots':
            ShowPlots = True
        elif opt == '--processes':
            nprocs = int(arg)

    # Check that the number of modes requested doesn't exceed 512
    if numModes > 512:
        print('Setting maximum number of expansion modes: 512.')
        numModes = 512

    # Check that only one configuration is chosen
    configs = [ExodusSingleConn, ExodusMultiConn, SCRIPwithoutConn, SCRIPwithConn]
    numConfigs = sum(bool(x) for x in configs)
    if numConfigs > 1:
        print('ONE mesh configuration option must be set!')
        print('None of the options are set.')
        sys.exit(2)

    if EvaluateAll:
        EvaluateTPW = EvaluateCFR = EvaluateTPO = EvaluateA1 = EvaluateA2 = True

    if 2*sampleOrder-1 < numModes:
        print("WARNING: The quadrature sampling order of %d is insufficient to exactly integrate SPH expansions of order %d!" % (
            sampleOrder, numModes))

    return sampleMesh, numModes, seed, \
        sampleCentroid, sampleOrder, \
        EvaluateTPW, EvaluateCFR, EvaluateTPO, \
        EvaluateA1, EvaluateA2, ShowPlots, \
        ExodusSingleConn, ExodusMultiConn, SCRIPwithoutConn, \
        SCRIPwithConn, SpectralElement, nprocs


if __name__ == '__main__':
    print('Welcome to CANGA remapping intercomparison field generator!')
    print('Authors: Jorge Guerra, Vijay Mahadevan, Paul Ullrich, 2019')

    # Parse the commandline! COMMENT OUT TO RUN IN IDE
    mesh_file, ND, seed, sampleCentroid, sampleOrder, \
        EvaluateTPW, EvaluateCFR, EvaluateTPO, \
        EvaluateA1, EvaluateA2, ShowPlots, \
        ExodusSingleConn, ExodusMultiConn, SCRIPwithoutConn, \
        SCRIPwithConn, SpectralElement, nprocs \
        = parseCommandLine(sys.argv[1:])

    # Set the name for the new data file
    stripDir = mesh_file.split('/')
    onlyFilename = stripDir[len(stripDir)-1]
    data_file = 'sample_NM' + str(ND) + '_O' + str(sampleOrder) + '_' + (onlyFilename.split('.'))[0]

    # Let us decipher what our final output file name should be with approrpriate suffixes
    outFileName = data_file

    if SpectralElement:
        outFileName += '_GLL'

    if EvaluateTPW:
        outFileName += '_TPW'
    if EvaluateCFR:
        outFileName += '_CFR'
    if EvaluateTPO:
        outFileName += '_TPO'
    if EvaluateA1:
        outFileName += '_A1'
    if EvaluateA2:
        outFileName += '_A2'
    outFileName += '.nc'

    print('File name for sampled mesh data: ', outFileName)
    print('Number of SH degrees for sampling set to: ', ND)
    print('Maximum Gaussian quadrature order to be used: ', 2*sampleOrder-1)

    if ExodusSingleConn or ExodusMultiConn:

        if SpectralElement:
            connCell = 'element_gll_conn'
            coordCell = 'grid_gll_cart'
        else:
            if ExodusSingleConn:
                connCell = 'connect1'
            elif ExodusMultiConn:
                connCell = 'connect0'

            coordCell = 'coord'

        # Open the .g mesh files for reading
        m_fid = Dataset(mesh_file)

        # Get connectivity and coordinate arrays (check for multiple connectivity)
        varCon = m_fid.variables[connCell][:]
        varCoord = m_fid.variables[coordCell][:]

        # Get the rectilinear attribute if available
        try:
            print('Rectilinear mesh detected; field variable written as 2D')
            rectilinear = m_fid.rectilinear
            # Get the 2D size of the field array from mesh file
            NLON = m_fid.rectilinear_dim1_size
            NLAT = m_fid.rectilinear_dim0_size
        except:
            print('NOT a rectilinear mesh.')
            rectilinear = False

    elif ExodusMultiConn:
        numElTypes = 'num_el_blk'
        numDims = 'cart_dims'
        connCell = 'element_corners_id'
        coordCell = 'grid_corners_cart'
        numVerts = 'grid_corners_size'

        # Open the .g mesh files for reading
        m_fid = Dataset(mesh_file)

        start = time.time()
        # Get connectivity and coordinate arrays
        varConnList = []
        numVertList = []
        numConnBlocks = len(m_fid.dimensions[numElTypes])
        for cc in range(numConnBlocks):
            # Get this connectivity array (El X corners)
            connName = 'connect' + str(cc+1)
            thisConn = m_fid.variables[connName][:]
            # Get the number of corners for this connectivity block
            numVertList.append(thisConn.shape[1])  # Column dimension of connectivity
            # Append to the list of connectivity blocks
            varConnList.append(m_fid.variables[connName][:])

        # Get the maximum number of vertices
        maxVerts = np.amax(np.array(numVertList))
        # Loop over the blocks again and pad columns up to the max vertices
        for cc in range(numConnBlocks):
            numVert2Pad = maxVerts - numVertList[cc]

            if numVert2Pad == 0:
                continue

            # Pad with redundant last coord ID up to the max vertices
            lastCol = np.expand_dims(varConnList[cc][:, -1], axis=1)
            thisPadding = np.matlib.repmat(lastCol, 1, numVert2Pad)
            varConnList[cc] = np.hstack((varConnList[cc], thisPadding))

        # Vertical stack of the connectivity lists
        varCon = np.vstack(tuple(varConnList))
        varCoord = m_fid.variables['coord'][:]

        try:
            print('Storing connectivity and coordinate arrays from Exodus mesh files.')
            numEdges = 'num_nod_per_el'
            numCells = 'num_el_in_blk'
            meshFileOut = m_fid.createDimension(numEdges, maxVerts)
            meshFileOut = m_fid.createDimension(numCells, varCon.shape[0])
            meshFileOut = m_fid.createDimension(numVerts, np.size(varCoord, 1))
            meshFileOut = m_fid.createDimension(numDims, 3)
            meshFileOut = m_fid.createVariable(connCell, 'i4', (numCells, numEdges))
            meshFileOut[:] = varCon
            meshFileOut = m_fid.createVariable(coordCell, 'f8', (numDims, numVerts))
            meshFileOut[:] = varCoord

        except RuntimeError:
            print('Cell connectivity and grid vertices exist in mesh data file.')

        endt = time.time()
        print('Time to precompute EXODUS multi-connectivity mesh info (sec): ', endt - start)

    elif SCRIPwithoutConn:
        numEdges = 'grid_corners'
        numCells = 'grid_size'
        numDims = 'cart_dims'
        numVerts = 'grid_corners_size'

        if SpectralElement:
            connCell = 'element_gll_conn'
            coordCell = 'grid_gll_cart'
        else:
            connCell = 'element_corners_id'
            coordCell = 'grid_corners_cart'

        # Open the .nc SCRIP files for reading
        m_fid = Dataset(mesh_file)

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

        if SpectralElement:
            connCell = 'element_gll_conn'
            coordCell = 'grid_gll_cart'
        else:
            connCell = 'element_corners'
            coordCell = 'grid_corners_cart'

        # Open the .nc SCRIP files for reading
        m_fid = Dataset(mesh_file)

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

    if SpectralElement:
        # Compute Lon/Lat coordinates from GLL nodes
        varLonLat = sphcrt.computeCart2LL(varCoord.T)
    else:
        # Compute Lon/Lat coordinates from centroids
        varCent = sphcrt.computeCentroids(varCon, varCoord)
        varLonLat = sphcrt.computeCart2LL(varCent)

    # Convert to degrees from radians
    varLonLat_deg = 180.0 / mt.pi * varLonLat

    m_fid.close()

    # Define our global variables for fields
    TPWvar = np.zeros(3)
    CFRvar = np.zeros(3)
    TPOvar = np.zeros(3)

    # %% Begin the SH reconstructions
    def Evaluate_TPW_Field():
        start = time.time()
        print('Computing Total Precipitable Water on sampling mesh...')
        # Set the power spectrum coefficients
        lfPower = [5.84729561e+04, -2.91678103e-04, -5.83966265e+04]
        hfPower = [2.17936330e+02, -1.99788552e+00, -7.94469251e-04]
        degIntersect = 1.8161917668847762
        # Compute the parent power spectrum for TPW
        degsTPW, psdTPW = computeSpectrum(ND, lfPower, hfPower, degIntersect)

        # Set the low degree coefficients (large scale structures)
        coeffsLD_TPW = np.array([[2.45709150e+01, 0.0, 0.0, 0.0],
                                 [4.00222122e+00, 2.39412571e+00, 0.0, 0.0],
                                 [-1.36433589e+01, 3.90520866e-03, 4.70350344e-01, 0.0],
                                 [-3.54931720e+00, -1.23629157e+00, 4.01454924e-01, 1.76782768e+00]])

        # Make the SH coefficients object for this field
        clmTPW = computeNormalizedCoefficients(ND, psdTPW, coeffsLD_TPW)

        # Evaluate actual spherical harmonic modes as solution;
        # change ls, ms below
        # lmax = 100
        # clmTPW = pyshtools.SHCoeffs.from_zeros(lmax)
        # clmTPW.set_coeffs(values=[1], ls=[2], ms=[2])

        # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
        # Expand the coefficients and check the field
        if sampleCentroid or SpectralElement:
            TPWvar = clmTPW.expand(lon=varLonLat_deg[:, 0], lat=varLonLat_deg[:, 1])
        else:
            TPWvar = computeCellAverage(clmTPW, varCon, varCoord, sampleOrder, True, nprocs)
            print('Total Precipitable Water Global integral: ', np.sum(TPWvar))

        # Compute rescaled data from 0.0 to max
        minTPW = np.amin(TPWvar)
        maxTPW = np.amax(TPWvar)
        deltaTPW = abs(maxTPW - minTPW)
        deltaTPW = deltaTPW if deltaTPW > 1e-10 else 1.0
        TPWvar = np.add(TPWvar, -minTPW)
        TPWvar *= maxTPW / deltaTPW
        endt = time.time()
        print('Time to compute TPW (mm): ', endt - start)

        return_dict['TPWvar'] = TPWvar

    # %%
    def Evaluate_CFR_Field():
        start = time.time()
        print('Computing Cloud Fraction on sampling mesh...')
        # Set the power spectrum coefficients
        lfPower = [8.38954430e+00, -1.85962382e-04, -8.38439294e+00]
        hfPower = [1.25594628e-01, -1.99203168e+00,  1.91763519e-06]
        degIntersect = 8.322269484619733
        # Compute the parent power spectrum for CFR
        degsCFR, psdCFR = computeSpectrum(ND, lfPower, hfPower, degIntersect)

        # Set the low degree coefficients (large scale structures)
        coeffsLD_CFR = np.array([[6.65795054e-01, 0.0, 0.0, 0.0],
                                 [-2.45480409e-02, 2.24697424e-02, 0.0, 0.0],
                                 [5.72322008e-02, 3.41184683e-02, -7.71082815e-03, 0.0],
                                 [1.86562455e-02, 4.34697733e-04, 8.91735978e-03, -5.53756958e-03]])

        # Make the SH coefficients object for this field
        clmCFR = computeNormalizedCoefficients(ND, psdCFR, coeffsLD_CFR)

        # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
        # Expand the coefficients and check the field
        if sampleCentroid or SpectralElement:
            CFRvar = clmCFR.expand(lon=varLonLat_deg[:, 0], lat=varLonLat_deg[:, 1])
        else:
            CFRvar = computeCellAverage(clmCFR, varCon, varCoord, sampleOrder, True, nprocs)
            print('Cloud Fraction Global integral: ', np.sum(CFRvar))

        # Compute rescaled data from 0.0 to max
        minCFR = np.amin(CFRvar)
        maxCFR = np.amax(CFRvar)
        deltaCFR = abs(maxCFR - minCFR)
        deltaCFR = deltaCFR if deltaCFR > 1e-10 else 1.0
        CFRvar = np.add(CFRvar, -minCFR)
        CFRvar *= maxCFR / deltaCFR
        #  Set all values greater than 1.0 to 1.0 (creates discontinuities)
        CFRvar[CFRvar >= 1.0] = 1.0

        endt = time.time()
        print('Time to compute CFR (0.0 to 1.0): ', endt - start)

        return_dict['CFRvar'] = CFRvar

    # %%
    def Evaluate_TPO_Field():
        start = time.time()
        print('Computing Global Terrain on sampling mesh...')
        # Set the power spectrum coefficients
        lfPower = [1.79242815e+05, -4.28193211e+01,  7.68040558e+05]
        hfPower = [9.56198160e+06, -1.85485966e+00, -2.63553217e+01]
        degIntersect = 3.8942282772035255
        # Compute the parent power spectrum for CFR
        degsTPO, psdTPO = computeSpectrum(ND, lfPower, hfPower, degIntersect)

        # Set the low degree coefficients (large scale structures)
        coeffsLD_TPO = np.array([[-2.38452711e+03, 0.0, 0.0, 0.0],
                                 [-6.47223253e+02, -6.06453097e+02, 0.0, 0.0],
                                 [5.67394318e+02, 3.32672611e+02, -4.17639577e+02, 0.0],
                                 [1.57403492e+02, 1.52896988e+02, 4.47106726e+02, -1.40553447e+02]])

        # Make the SH coefficients object for this field
        clmTPO = computeNormalizedCoefficients(ND, psdTPO, coeffsLD_TPO)

        # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
        # Expand the coefficients and check the field
        if sampleCentroid or SpectralElement:
            TPOvar = clmTPO.expand(lon=varLonLat_deg[:, 0], lat=varLonLat_deg[:, 1])
        else:
            TPOvar = computeCellAverage(clmTPO, varCon, varCoord, sampleOrder, True, nprocs)
            print('Global Terrain Global integral: ', np.sum(TPOvar))

        # Rescale to -1.0 to 1.0
        minTPO = np.amin(TPOvar)
        maxTPO = np.amax(TPOvar)
        deltaTPO = abs(maxTPO - minTPO)
        deltaTPO = deltaTPO if deltaTPO > 1e-10 else 1.0
        TPOvar = np.add(TPOvar, -0.5 * (maxTPO + minTPO))
        TPOvar *= 2.0 / deltaTPO

        # Rescale topography to real Earth max/min
        minTPO = -10994.0  # Depth at Challenger Deep
        maxTPO = 8848.0  # Elevation of Mt. Everest ASL
        deltaTPO = abs(maxTPO - minTPO)
        TPOvar *= (0.5 * deltaTPO)
        TPOvar += 0.5 * (maxTPO + minTPO)

        endt = time.time()
        print('Time to compute TPO (m): ', endt - start)

        return_dict['TPOvar'] = TPOvar

    # %%
    def Evaluate_A1_Field():
        start = time.time()
        print('Computing Analytical Field 1 sampling on mesh...')

        # Evaluate actual spherical harmonic modes as solution;
        # change ls, ms below
        lmax = 100
        clmA1 = pyshtools.SHCoeffs.from_zeros(lmax)
        # This evaluates P_3^3
        clmA1.set_coeffs(values=[1], ls=[3], ms=[2])
        clmA1.set_coeffs(values=[1], ls=[3], ms=[3])

        # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
        # Expand the coefficients and check the field
        if sampleCentroid or SpectralElement:
            A1var = clmA1.expand(lon=varLonLat_deg[:, 0], lat=varLonLat_deg[:, 1])
            print('Analytical Solution 1 Global sum: ', np.sum(A1var)/A1var.shape[0])
        else:
            A1var = computeCellAverage(clmA1, varCon, varCoord, sampleOrder, True, nprocs)
            print('Analytical Solution 1 Global integral: ', np.sum(A1var)/A1var.shape[0])

        endt = time.time()
        print('Time to compute A1 Field: ', endt - start)

        return_dict['A1var'] = A1var

    # %%
    def Evaluate_A2_Field():
        start = time.time()
        print('Computing Analytical Field 2 sampling on mesh...')

        # THIS NEEDS TO CHANGE TO SUPPORT FE GRIDS
        # Expand the coefficients and check the field
        # if sampleCentroid or SpectralElement:
        if sampleCentroid or SpectralElement:
            A2var = evaluate_field_a2(lon=varLonLat[:, 0], lat=varLonLat[:, 1])
            print('Analytical Solution 2 Global sum: ', np.sum(A2var)/A2var.shape[0])
        else:
            # A2var = computeCellAverageSerial(evaluate_field_a2, varCon, varCoord, sampleOrder, True)
            A2var = computeCellAverage(evaluate_field_a2, varCon, varCoord, sampleOrder, True, nprocs)
            print('Analytical Solution 2 Global integral: ', np.sum(A2var)/A2var.shape[0])

        endt = time.time()
        print('Time to compute A2 Field: ', endt - start)

        return_dict['A2var'] = A2var

    # %%

    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    # Let us aggregate all the jobs that need to be done and then
    # let the multiprocessing manager take care of it.
    jobs = []
    evaluation_routines = []
    if EvaluateTPW:
        evaluation_routines.append(Evaluate_TPW_Field)
    if EvaluateCFR:
        evaluation_routines.append(Evaluate_CFR_Field)
    if EvaluateTPO:
        evaluation_routines.append(Evaluate_TPO_Field)
    if EvaluateA1:
        evaluation_routines.append(Evaluate_A1_Field)
    if EvaluateA2:
        evaluation_routines.append(Evaluate_A2_Field)
    for fn in evaluation_routines:
        p = Process(target=fn)
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()

    if EvaluateTPW:
        TPWvar = return_dict['TPWvar']
    if EvaluateCFR:
        CFRvar = return_dict['CFRvar']
    if EvaluateTPO:
        TPOvar = return_dict['TPOvar']
    if EvaluateA1:
        A1var = return_dict['A1var']
    if EvaluateA2:
        A2var = return_dict['A2var']

    # %% Copy grid files and store the new test data (source and target)
    shutil.copy(mesh_file, outFileName)

    # write lon, lat, and test data variables
    data_fid = Dataset(outFileName, 'a')

    # Set the dimension name depending on the mesh file format
    if ExodusSingleConn:
        numCells = 'num_el_in_blk1'
    elif ExodusMultiConn:
        numCells = 'num_el_in_blk0'
    elif SCRIPwithoutConn:
        numCells = 'grid_size'
    elif SCRIPwithConn:
        numCells = 'ncells'

    if SpectralElement:
        numCells = 'grid_gll_size'

    # Process the sampling file
    if SCRIPwithConn:
        lonNC = data_fid.createVariable('nlon', 'f8', (numCells,))
        lonNC[:] = varLonLat_deg[:, 0]
        latNC = data_fid.createVariable('nlat', 'f8', (numCells,))
        latNC[:] = varLonLat_deg[:, 1]
    else:
        lonNC = data_fid.createVariable(
            'lon', 'f8', (numCells,)) if 'lon' not in data_fid.variables.keys() else data_fid.variables['lon']
        lonNC[:] = varLonLat_deg[:, 0]
        latNC = data_fid.createVariable(
            'lat', 'f8', (numCells,)) if 'lat' not in data_fid.variables.keys() else data_fid.variables['lat']
        latNC[:] = varLonLat_deg[:, 1]

    if rectilinear:
        slon = 'lonDim'
        slat = 'latDim'
        data_fid.createDimension(slon, NLON)
        data_fid.createDimension(slat, NLAT)

        if EvaluateTPW:
            TPWNC = data_fid.createVariable('TotalPrecipWater', 'f8', (slat, slon)) if 'TotalPrecipWater' not in data_fid.variables.keys(
            ) else data_fid.variables['TotalPrecipWater']
            field = np.reshape(TPWvar, (NLAT, NLON))
            TPWNC[:] = field
        if EvaluateCFR:
            CFRNC = data_fid.createVariable('CloudFraction', 'f8', (slat, slon)) if 'CloudFraction' not in data_fid.variables.keys(
            ) else data_fid.variables['CloudFraction']
            field = np.reshape(CFRvar, (NLAT, NLON))
            CFRNC[:] = field
        if EvaluateTPO:
            TPONC = data_fid.createVariable(
                'Topography', 'f8', (slat, slon)) if 'Topography' not in data_fid.variables.keys() else data_fid.variables['Topography']
            field = np.reshape(TPOvar, (NLAT, NLON))
            TPONC[:] = field
        if EvaluateA1:
            A1NC = data_fid.createVariable('AnalyticalFun1', 'f8', (slat, slon)) if 'AnalyticalFun1' not in data_fid.variables.keys(
            ) else data_fid.variables['AnalyticalFun1']
            field = np.reshape(A1var, (NLAT, NLON))
            A1NC[:] = field
        if EvaluateA2:
            A2NC = data_fid.createVariable('AnalyticalFun2', 'f8', (slat, slon)) if 'AnalyticalFun2' not in data_fid.variables.keys(
            ) else data_fid.variables['AnalyticalFun2']
            field = np.reshape(A2var, (NLAT, NLON))
            A2NC[:] = field
    else:
        if EvaluateTPW:
            TPWNC = data_fid.createVariable('TotalPrecipWater', 'f8', (numCells,)) if 'TotalPrecipWater' not in data_fid.variables.keys(
            ) else data_fid.variables['TotalPrecipWater']
            TPWNC[:] = TPWvar
        if EvaluateCFR:
            CFRNC = data_fid.createVariable('CloudFraction', 'f8', (numCells,)) if 'CloudFraction' not in data_fid.variables.keys(
            ) else data_fid.variables['CloudFraction']
            CFRNC[:] = CFRvar
        if EvaluateTPO:
            TPONC = data_fid.createVariable(
                'Topography', 'f8', (numCells,)) if 'Topography' not in data_fid.variables.keys() else data_fid.variables['Topography']
            TPONC[:] = TPOvar
        if EvaluateA1:
            A1NC = data_fid.createVariable('AnalyticalFun1', 'f8', (numCells,)) if 'AnalyticalFun1' not in data_fid.variables.keys(
            ) else data_fid.variables['AnalyticalFun1']
            A1NC[:] = A1var
        if EvaluateA2:
            A2NC = data_fid.createVariable('AnalyticalFun2', 'f8', (numCells,)) if 'AnalyticalFun2' not in data_fid.variables.keys(
            ) else data_fid.variables['AnalyticalFun2']
            A2NC[:] = A2var

    # Close the files out.
    data_fid.close()

    # '''
    # %% Check the data with triangular surface plot
    if ShowPlots:
        points2D = varLonLat
        tri = Delaunay(points2D)
        simplices = tri.simplices

        # %% Plot Total Precipitable Water
        if EvaluateTPW:
            fig1 = FF.create_trisurf(x=varLonLat[:, 0], y=varLonLat[:, 1], z=TPWvar, height=800, width=1200,
                                     simplices=simplices, colormap="Portland", plot_edges=False,
                                     title="Total Precipitable Water Check (mm)", aspectratio=dict(x=1, y=1, z=0.3))
            py.offline.plot(fig1, filename='TPW' + data_file + '.html')

            '''
                     import matplotlib.pyplot as plt
                     from mpl_toolkits.basemap import Basemap

                     fig = plt.figure(figsize=(10, 8))
                     m = Basemap(projection='lcc', resolution='c',
                            width=8E6, height=8E6,
                            lat_0=45, lon_0=-100,)
                     m.shadedrelief(scale=0.5)
                     m.pcolormesh(varLonLat[:,0], varLonLat[:,1], TPWvar,
                            latlon=True, cmap='RdBu_r')
                     plt.clim(-8, 8)
                     m.drawcoastlines(color='lightgray')

                     plt.title('January 2014 Temperature Anomaly')
                     plt.colorbar(label='temperature anomaly (°C)');
                     '''

        # %% Plot Cloud Fraction
        if EvaluateCFR:
            fig1 = FF.create_trisurf(x=varLonLat[:, 0], y=varLonLat[:, 1], z=CFRvar, height=800, width=1200,
                                     simplices=simplices, colormap="Portland", plot_edges=False,
                                     title="Cloud Fraction Check (0.0-1.0)", aspectratio=dict(x=1, y=1, z=0.3))
            py.offline.plot(fig1, filename='CFR' + data_file + '.html')
        # %% Plot Topography
        if EvaluateTPO:
            fig1 = FF.create_trisurf(x=varLonLat[:, 0], y=varLonLat[:, 1], z=TPOvar, height=800, width=1200,
                                     simplices=simplices, colormap="Portland", plot_edges=False,
                                     title="Global Topography (m)", aspectratio=dict(x=1, y=1, z=0.3))
            py.offline.plot(fig1, filename='TPO' + data_file + '.html')
        # '''
        # %% Plot Topography
        if EvaluateA1:
            fig1 = FF.create_trisurf(x=varLonLat[:, 0], y=varLonLat[:, 1], z=A1var, height=800, width=1200,
                                     simplices=simplices, colormap="Portland", plot_edges=False,
                                     title="Analytical Function 1 (SPH(3,3))", aspectratio=dict(x=1, y=1, z=0.3))
            py.offline.plot(fig1, filename='A1' + data_file + '.html')
        # '''
        # %% Plot Topography
        if EvaluateA2:
            fig1 = FF.create_trisurf(
                x=varLonLat[:, 0],
                y=varLonLat[:, 1],
                z=A2var, height=800, width=1200, simplices=simplices, colormap="Portland", plot_edges=False,
                title="Analytical Function 2: (2.0 + cos^2(lat) * cos(2.0 * lon))",
                aspectratio=dict(x=1, y=1, z=0.3))
            py.offline.plot(fig1, filename='A2' + data_file + '.html')
        # '''

        # %% Check the evaluated spectra
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


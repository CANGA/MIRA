#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  8 10:16:17 2019

 Creates (reconstructs) a coordinate array and connectivity based on SCRIP lon/lat data
 Generates grid ID's on a per cell basis without reordering. This function has a hard
 coded tolerance for coincident nodes. USE/SET WITH CARE.

@author: TempestGuerra
"""

import numpy as np
from scipy.spatial import cKDTree


def computeCoordConFastSCRIP(lon, lat):

    # Initialize and set tolerance for coincident grids
    kdleafs = 64
    COINCIDENT_TOLERANCE = 1.0E-14
    NC = np.size(lon, axis=0)
    NG = np.size(lon, axis=1)
    # Coordinate array starts as two triplets (REMOVED AT THE END)
    gridCoord = np.zeros((1, 6))
    # Connectivity array same size as the input
    cellCon = np.zeros((NC, NG))
    # Initialize gridID
    gridID = 0
    # Loop over the raw connectivity cells
    for ii in range(NC):

        # Loop over grids in the connectivity
        for jj in range(NG):

            # Set each grid as a new grid
            gridID += 1

            # Get a grid at the current cell (unit radius)
            thisGrid = np.array(
                [gridID, ii + 1, jj + 1, lon[ii, jj], lat[ii, jj], 1.0])

            # Put the new grid into the raw coordinate array
            gridCoord = np.append(gridCoord, [thisGrid.T], axis=0)

            # Put this grid ID into the connectivity
            cellCon[ii, jj] = int(gridID)

    # Remove the first triplet from the coordinate array
    gridCoord = np.delete(gridCoord, 0, axis=0)

    # Make a copy of the gridCoord to relabel coincidents
    newGridCoord = gridCoord

    # Build the KDtree for the coordinate array (using lon/lat only)
    pointList = gridCoord[:, 3:5]
    # Compute the KD tree object
    pointTree = cKDTree(pointList, leafsize=kdleafs)

    # Find the coincident nodes for each query node
    NV = np.size(gridCoord, axis=0)
    for ii in range(NV):
        # Get this gridID
        gridID = gridCoord[ii, 0]
        # Compute the list of nodes that are coincident
        thisGrid = pointList[ii, :]
        ndex = pointTree.query_ball_point(
            thisGrid, COINCIDENT_TOLERANCE, p=2, eps=0)

        NN = len(ndex)
        # Loop over the coincidents
        for cc in range(NN):
            # Overwrite new grid array with gridID
            newGridCoord[ndex[cc], 0] = gridID

    # Sort the new grid coordinate array
    sortDex = np.argsort(newGridCoord[:, 0])
    sortedNewGridCoord = newGridCoord[sortDex, :]

    # Renumber column of gridIDs to be in ascending order
    kk = 1
    for ii in range(1, NV):
        # Check for coincidence with previous node
        if sortedNewGridCoord[ii, 0] == sortedNewGridCoord[ii - 1, 0]:
            sortedNewGridCoord[ii, 0] = sortedNewGridCoord[ii - 1, 0]
        else:
            kk += 1
            sortedNewGridCoord[ii, 0] = kk

        # Fix the connectivity with sorted and renumbered nodes
        rdex = int(sortedNewGridCoord[ii, 1] - 1)
        cdex = int(sortedNewGridCoord[ii, 2] - 1)
        cellCon[rdex, cdex] = sortedNewGridCoord[ii, 0]

    # Trim the grid coordinate array of coincident nodes
    cleanNewGridCoord = np.zeros((1, 4))
    # This index slice selects for [ID lon lat 1.0]
    adex = [0, 3, 4, 5]
    cleanNewGridCoord[0, :] = sortedNewGridCoord[0, adex]
    for ii in range(1, NV):
        # Check for coincidence with previous node
        if sortedNewGridCoord[ii, 0] != sortedNewGridCoord[ii - 1, 0]:
            thisGrid = sortedNewGridCoord[ii, adex]
            cleanNewGridCoord = np.append(
                cleanNewGridCoord, [thisGrid], axis=0)

    return cleanNewGridCoord, cellCon.astype(int)

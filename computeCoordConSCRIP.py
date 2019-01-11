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

def computeKDtreeSearch(thisGrid, pointList, COINCIDENT_TOLERANCE):
       
       newPointList = pointList
       # compute the KD tree object
       pointTree = cKDTree(pointList, leafsize=20)
       # compute a few (10) of the nearest neighbors 
       NN = 10
       dd, ndex = pointTree.query(thisGrid, k=range(1,NN), eps=0, p=2, \
                                distance_upper_bound=np.Inf, n_jobs=1)
       
       # Get the minimum distance from the set of nearest neighbors
       dmin = np.amin(dd)
       
       # If a coincident node is among the nearest neighbors, set the gridID
       if dmin <= COINCIDENT_TOLERANCE:
              gdex = np.argmin(dd)
              gridID = ndex[gdex]
       else:
              # This is a new grid to the coordinate list, set the gridID
              newPointList = np.append(pointList, [thisGrid.T], axis=0)
              gridID = np.size(newPointList, axis=0)
              
       
       return gridID, newPointList

def computeLinearSearch(thisGrid, pointList, NS, NE, COINCIDENT_TOLERANCE):
       gridID = 0
       newGrid = False
       newPointList = pointList
              
       # Check thisGrid to see if it is in the point list
       for kk in range(NS, NE, 1):
              
              pDiff = np.subtract(thisGrid.T, pointList[kk,:])
              
              # Get out of this search if coincidence is met
              if np.sum(np.power(pDiff, 2)) <= COINCIDENT_TOLERANCE:
                     gridID = kk
                     newGrid = False
                     break;
              else:
                     newGrid = True
                     gridID = kk + 1
       
       if newGrid:
              newPointList = np.append(pointList, [thisGrid.T], axis=0)
       
       return gridID, newPointList
       

def computeCoordConSCRIP(lon, lat):
       
       # Initialize and set tolerance for coincident grids
       INITIAL_SEARCH = 100
       COINCIDENT_TOLERANCE = 1.0E-14
       NC = np.size(lon, axis=0)
       NG = np.size(lon, axis=1)
       # Coordinate array starts as two triplets (REMOVED AT THE END) 
       varCoord = np.zeros((1,3))
       # Connectivity array same size as the input
       varCon = np.zeros((NC, NG))
       
       # Loop over the raw connectivity cells
       for ii in range(NC):
              
              # Loop over grids in the connectivity
              for jj in range(NG):
                     
                     # Get a grid at the current cell (unit radius)
                     thisGrid = np.array([lon[ii,jj], lat[ii,jj], 1.0])
                     
                     # Linear search on the first few cells
                     if ii < INITIAL_SEARCH:
                            NP = np.size(varCoord, axis=0)
                            gridID, varCoord = computeLinearSearch(thisGrid, varCoord, 0, NP, COINCIDENT_TOLERANCE)
                     else:
                            # Use a KD tree in each search
                            gridID, varCoord = computeKDtreeSearch(thisGrid, varCoord, COINCIDENT_TOLERANCE)
                            
                     # Put this grid ID into the connectivity
                     varCon[ii,jj] = int(gridID)
                     
              #print(ii, varCon[ii,:])
              
       # Remove the first triplet from the coordinate array
       varCoord = np.delete(varCoord, 0, axis=0)
       
       return varCoord, varCon
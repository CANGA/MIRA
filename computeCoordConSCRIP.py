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
import random as rand

def computeSpacePartitionSearch(thisGrid, pointList, COINCIDENT_TOLERANCE):
       
       gridID = 0
       NP = np.size(pointList, axis=0)
       newPointList = pointList
              
       # Compute a partion into 4... 
       aloc = int(NP / 4)
       bloc = int(NP / 2)
       cloc = int(3 * NP / 4)
       # Pick random indices in each partition
       pdex1 = rand.randint(0, aloc)
       pdex2 = rand.randint(aloc, bloc)
       pdex3 = rand.randint(bloc, cloc)
       pdex4 = rand.randint(cloc, NP - 1)
       # Compute differences
       diff1 = np.subtract(thisGrid.T, pointList[pdex1,:])
       diff2 = np.subtract(thisGrid.T, pointList[pdex2,:])
       diff3 = np.subtract(thisGrid.T, pointList[pdex3,:])
       diff4 = np.subtract(thisGrid.T, pointList[pdex4,:])
       # Compute distances
       dist1 = np.linalg.norm(diff1)
       dist2 = np.linalg.norm(diff2)
       dist3 = np.linalg.norm(diff3)
       dist4 = np.linalg.norm(diff4)
       dists = np.array([dist1, dist2, dist3, dist4])
       
       # Compute minimum of branch distances and get the winning partition
       mdex = 1
       
       # Compare distances and search only in half the points
       if mdex == 1:
              NS = 0
              NE = int(NP / 2)
       elif dist2 < dist1:
              NS = int(NP / 2)
              NE = NP
       else:
              NS = 0
              NE = NP
              
       # Compute linear search over the winning partition
       gridID, newPointList = computeLinearSearch(thisGrid, pointList, NS, NE, COINCIDENT_TOLERANCE)
       
       return gridID, newPointList

def computeLinearSearch(thisGrid, pointList, NS, NE, COINCIDENT_TOLERANCE):
       gridID = 0
       newGrid = False
       newPointList = pointList
              
       # Check thisGrid to see if it is in the point list
       for kk in range(NS, NE, 1):
              
              pDiff = np.subtract(thisGrid.T, pointList[kk,:])
              
              # Get out of this search if coincidence is met
              if np.linalg.norm(pDiff) <= COINCIDENT_TOLERANCE:
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
       COINCIDENT_TOLERANCE = 1.0E-12
       NC = np.size(lon, axis=0)
       NG = np.size(lon, axis=1)
       # Coordinate array starts as two triplets (REMOVED AT THE END) 
       varCoord = np.zeros((1,3))
       # Connectivity array same size as the input
       varCon = np.zeros((NC, NG))
       
       # Loop over the raw connectivity cells
       NC = 1000
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
                            # Simple partitioned search for coincident node
                            gridID, varCoord = computeSpacePartitionSearch(thisGrid, varCoord, COINCIDENT_TOLERANCE)
                            
                     # Put this grid ID into the connectivity
                     varCon[ii,jj] = gridID
       
       return varCoord, varCon
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:26:47 2018

Compute a low order FV gradient on the manifold based on the adjacency stencil 
around a cell. Corresponds to "strategy 3" in Barth & Jespersen, 1989 page 6.

By construction, all the vectors involved start at the origin... convenient.
Parameterized flux integral assumes constant radius and field value.

@author: jeguerra
"""

import numpy as np
import math as mt
from computeAreaIntegral import computeAreaIntegral
import computeSphericalCartesianTransforms as sphcrt

import multiprocessing
from multiprocessing import Process
from itertools import repeat


def computeGradientFV3_Private(NC, NP, jj, varField, varCon, varCoords, varStenDex, radius, cellCoords):
      SF = np.float64
       
      # Gradients are 3 component vectors
      nc = 3
       
      # Set the main return variable
      # varGradient = np.zeros(nc, dtype=SF)

      # Loop over the cells
      # for jj in range(NC):
      pdex = np.array(range(NP), dtype = int)
      
      # Check for local degeneracy in stencil and fix connectivity
      for pp in range(NP):
              # Look for -1 in the adjacency stencil
              if varStenDex[jj,pp] <= 0:
                    pdex = np.delete(pdex, pp)
              else:
                    continue
              
      # Make pdex periodic
      pdexp = np.append(pdex, pdex[0])
      
      # Fetch the modified stencil
      thisStencil = varStenDex[jj,pdexp]
      
      # Initialize the new convex hull stencil
      convHullSten = varStenDex[jj,pdex]
      
      # Build the convex hull of cells around this cell
      for pp in range(len(convHullSten)):
              # Fetch consecutive pairs of cell id
              cid1 = thisStencil[pp] - 1
              cid2 = thisStencil[pp+1] - 1
              
              # Fetch consecutive pairs of stencils
              stn1 = varStenDex[cid1.astype(int),:]
              stn2 = varStenDex[cid2.astype(int),:]
              
              # Get the set intersection
              commonIds = list(set(stn1).intersection(stn2))
              # Get the common cell that is NOT the current target cell
              newCellId = [x for x in commonIds if x != jj+1]
              
              # Check new cell ID to be of length 1
              if len(newCellId) != 1:
                    # print('Found no neighboring cell or multiples in stencil!')
                    # print('New cell will NOT be recorded to the convex hull at cell: ', jj+1)
                    continue
              
              # Insert the new cell ID
              np.insert(convHullSten, pp+1, newCellId[0].astype(int))
      
      # Loop over the convex hull stencil and get dual edges map
      NS = len(convHullSten)
      fluxIntegral = np.zeros(nc, dtype=SF)
      dualEdgeMap = np.zeros((nc,NS))
      boundaryNorm = np.zeros((nc,NS))
      boundaryAngles = np.zeros((NS,1))
      for pp in range(NS):
              # Fetch the dual edge and store
              sid1 = convHullSten[pp] - 1
              sid1 = sid1.astype(int)
              # Make the dual polygon convex
              if pp == len(pdex) - 1:
                    sid2 = convHullSten[0] - 1
              else:
                    sid2 = convHullSten[pp+1] - 1
              sid2 = sid2.astype(int)
              
              # Store the dual mesh polygon
              dualEdgeMap[:,pp] = cellCoords[:,sid1]
              
              # Compute angles spanned by each boundary segment
              RE = 0.5 * (radius[sid1] + radius[sid2])
              unCoord1 = 1.0 / radius[sid1] * cellCoords[:,sid1]
              unCoord2 = 1.0 / radius[sid2] * cellCoords[:,sid2]
              boundaryAngles[pp] = mt.acos(np.dot(unCoord1, unCoord2))
              boundaryAngles[pp] = abs(boundaryAngles[pp])
              
              # Compute the stencil boundary normals
              boundaryNorm[:,pp] = np.cross(cellCoords[:,sid2], \
                                            cellCoords[:,sid1])
              bnMag = np.linalg.norm(boundaryNorm[:,pp])
              boundaryNorm[:,pp] = 1.0 / bnMag * boundaryNorm[:,pp]
              
              # Compute the weighted average of the two cell values AT the shared edge location
              vWeight = 0.5 * boundaryAngles[pp] * RE
              varAvg = varField[sid1] + varField[sid2]
              
              # Compute the integral over this edge
              fluxIntegral = np.add(fluxIntegral, \
                                vWeight * varAvg * \
                                boundaryNorm[:,pp])
                
      # Compute the dual polygon area
      areaD = computeAreaIntegral(None, dualEdgeMap, 6, False, True)
      
      # Compute the local gradient at this cell
      varGradient = 1.0 / areaD * fluxIntegral

      return varGradient#, areaD

def computeGradientFV3(varField, varCon, varCoords, varStenDex):
    SF = np.float64
    
    # Gradients are 3 component vectors
    nc = 3
    
    # Set the main return variable
    
    cellCoords = np.zeros((nc, varCon.shape[0]), dtype=SF)
    
    NC = int(varStenDex.shape[0])
    NP = int(varStenDex.shape[1])
    # areaD = np.zeros(NC, dtype=SF)
    
    # Precompute the cell centroid map
    cellCoords = sphcrt.computeCentroids(varCon, varCoords)
    radius = np.linalg.norm(cellCoords, axis=1)
    
    # Loop over each cell and get cell average
    pool = multiprocessing.Pool(processes=64)
    results = pool.starmap(computeGradientFV3_Private, zip(repeat(NC), repeat(NP), range(NC), repeat(varField), repeat(varCon), repeat(varCoords), repeat(varStenDex), repeat(radius), repeat(cellCoords)))
    pool.close()
    pool.join()
    #print(results)
    varGradient = np.reshape(np.array(results, dtype='f8').T, (nc, varField.shape[0]))
    #areaD  = np.array(results, dtype='f8')[:, 1]

    # print(varGradient)

    return varGradient



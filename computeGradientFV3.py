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
from computeCentroid import computeCentroid

import multiprocessing
from multiprocessing import Process
from itertools import repeat
from numba import jit
import time

@jit(nopython=True,parallel=False)
def loop_intersection(lst1, lst2):
    result = []
    for element1 in lst1:
        for element2 in lst2:
            if element1 == element2:
                result.append(element1)
    return result

@jit(nopython=True,parallel=False)
def gradientLoop(varGradient, fluxIntegralTotal, varField, areaInvD, NS, sid, fluxIntegral):
       NC = len(areaInvD)
       for elem in range(NC):
           optimizedComputeGradientFV3_Private(varGradient[:,elem], fluxIntegralTotal, varField, areaInvD[elem], sid[elem,0:NS[elem],:], fluxIntegral[elem,:,0:NS[elem]])

@jit(nopython=True,parallel=False)
def optimizedComputeGradientFV3_Private(varGradient, fluxIntegralTotal, varField, areaInvD, sid, fluxIntegral):

      # Gradients are 3 component vectors
      nc = 3

      # Loop over the convex hull stencil and get dual edges map
      NS = len(sid[:,0])

      # set entries to zero since this vector is reused by each call in the loop
      for i in range(nc):
          fluxIntegralTotal[i]=0
      for pp in range(NS):

            # Compute the weighted average of the two cell values AT the shared edge location
            # VSM: we have precomputed and stored values of sid already
            varAvg = varField[sid[pp,0]] + varField[sid[pp,1]]

            # Compute the integral over this edge as a weighted addition with precomputed
            # fluxIntegral on this edge
            fluxIntegralTotal = np.add(fluxIntegralTotal, varAvg * fluxIntegral[:,pp])

      # Compute the local gradient at this cell
      for i in range(nc):
          varGradient[i] = areaInvD * fluxIntegralTotal[i]

@jit(nopython=True,parallel=False)
def buildConvexHullOfCellsAroundCell(jj, NP, pdex, pdexp, varStenDex, convHullSten, thisStencil):
      # Check for local degeneracy in stencil and fix connectivity
      count = 0 
      for pp in range(NP):
            # Look for -1 in the adjacency stencil
            if varStenDex[jj,pp] <= 0:
                  count += 1
            else:
                  continue
      pdex_size = NP-count
      count = 0
      for pp in range(NP):
            # Look for -1 in the adjacency stencil
            if varStenDex[jj,pp] <= 0:
                  continue
            else:
                  pdex[count] = pp
                  count += 1
      # Make pdex periodic
      pdexp[pdex_size] = pdex[0]

      # Fetch the modified stencil
      for pp in range(pdex_size+1):
            thisStencil[pp] = varStenDex[jj,pdexp[pp]]
      for pp in range(pdex_size):
            convHullSten[pp] = varStenDex[jj,pdex[pp]]

      ## Build the convex hull of cells around this cell
      inserted_so_far = 0
      for pp in range(pdex_size):
            #print(NP, pdex_size)
            # Fetch consecutive pairs of cell id
            cid1 = thisStencil[pp] - 1
            cid2 = thisStencil[pp+1] - 1

            # Fetch consecutive pairs of stencils
            stn1 = varStenDex[cid1,:]
            stn2 = varStenDex[cid2,:]

            # Get the set intersection
            commonIds = loop_intersection(stn1, stn2)
            newCellId = [x for x in commonIds if x != jj+1]

            # Check new cell ID to be of length 1
            if len(newCellId) != 1:
                  # print('Found no neighboring cell or multiples in stencil!')
                  # print('New cell will NOT be recorded to the convex hull at cell: ', jj+1)
                  continue

            # Insert the new cell ID
            convHullSten[pdex_size+inserted_so_far] = newCellId[0]
            inserted_so_far+=1

      return (pdex_size, pdex_size+inserted_so_far)

@jit(nopython=True,parallel=False)
def computeDualEdgesMap(NS, convHullSten, pdex, dualEdgeMap, radius, cellCoords, boundaryAngles, boundaryNorm, sid, fluxIntegral):
      for pp in range(NS):
            # Fetch the dual edge and store
            sid1 = convHullSten[pp] - 1
            # Make the dual polygon convex
            if pp == len(pdex) - 1:
                  sid2 = convHullSten[0] - 1
            else:
                  sid2 = convHullSten[pp+1] - 1

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

            sid[pp, 0] = sid1
            sid[pp, 1] = sid2

            fluxIntegral[:,pp] =  vWeight * boundaryNorm[:,pp]

def precomputeGradientFV3Data_Private(NC, NP, jj, varCon, varCoords, varStenDex, radius, cellCoords):
      SF = np.float64

      # Gradients are 3 component vectors
      nc = 3

      pdex = np.zeros(shape=(NP,), dtype = int)
      pdexp = np.zeros(shape=(NP+1,), dtype = int)
      convHullSten = np.zeros(shape=(2*NP,), dtype = int)
      thisStencil = np.zeros(shape=(NP+1,), dtype = int)

      pdex_size, convHullSten_size = buildConvexHullOfCellsAroundCell(jj, NP, pdex, pdexp, varStenDex, convHullSten, thisStencil)

      # shrink from max allocation based on pdex_size
      pdex = pdex[:pdex_size]
      pdexp = pdexp[:pdex_size+1]
      thisStencil = thisStencil[:pdex_size+1]
      convHullSten = convHullSten[:convHullSten_size]

      NS = len(convHullSten)
      dualEdgeMap = np.zeros((nc,NS))
      boundaryNorm = np.zeros((nc,NS))
      boundaryAngles = np.zeros((NS,))
      sid = np.zeros((NS,2))
      fluxIntegral = np.zeros((nc,NS), dtype=SF)

      # Loop over the convex hull stencil and get dual edges map
      computeDualEdgesMap(NS, convHullSten, pdex, dualEdgeMap, radius, cellCoords, boundaryAngles, boundaryNorm, sid, fluxIntegral)

      # Compute the dual polygon area
      areaD = computeAreaIntegral(None, dualEdgeMap, 6, False, True)

      #############
      #  VSM: 
      # Store 1/area, [sid1, sid2], fluxIntegral[pp]
      #   where fluxIntegral[pp] = [vWeight * boundaryNorm[:,pp]]]
      # Then can compute: varAvg[pp] = varField[sid1] + varField[sid2]
      # And use the weighted addition:  fluxIntegralTotal = np.dot(fluxIntegral, varAvg) 
      # to compute the gradient as gradient = fluxIntegralTotal / area
      #############
      return 1.0 / areaD, sid, fluxIntegral

class ComputeGradientFV:

      def __init__(self, ctx, varCon, varCoords, varStenDex, NPROCS=4):
            ## Precompute key datastructures
            self.varCon = varCon
            self.varCoords = varCoords
            self.varStenDex = varStenDex
            self.NPROCS = NPROCS
            self.context = ctx

            # Precompute the cell centroid map
            self.cellCoords = sphcrt.computeCentroids(varCon, varCoords).T
            self.radius = np.linalg.norm(self.cellCoords, axis=0)
            NC = int(self.varStenDex.shape[0])
            self.areaInvD = np.zeros(NC, dtype=np.float64)
            self.cachedData = False

      def precomputeGradientFV3Data(self):

            if not self.cachedData:
                  SF = np.float64
                  NC = int(self.varStenDex.shape[0])
                  NP = int(self.varStenDex.shape[1])

                  pool = multiprocessing.Pool(processes=self.NPROCS)
                  results = pool.starmap(precomputeGradientFV3Data_Private, zip(repeat(NC), repeat(NP), range(NC), repeat(self.varCon), repeat(self.varCoords), repeat(self.varStenDex), repeat(self.radius), repeat(self.cellCoords)))
                  pool.close()
                  pool.join()

                  ## Store the data in the local object for future use
                  self.areaInvD = np.zeros(NC, dtype=SF)
                  self.NS = np.zeros(shape=(NC,), dtype='i4')
                  for elem in range(NC):
                        self.areaInvD[elem] = results[elem][0]
                        self.NS[elem] = np.array(results[elem][1], dtype='i4').shape[0]

                  max_NS = np.amax(self.NS)
                  nc = 3
                  self.sid = np.zeros(shape=(NC,max_NS,2), dtype='i4')
                  self.fluxIntegral = np.zeros((NC,nc,max_NS), dtype=np.float64)
                  for elem in range(NC):
                        self.sid[elem,0:self.NS[elem],:] = np.array(results[elem][1], dtype='i4')
                        self.fluxIntegral[elem,:,0:self.NS[elem]] = results[elem][2]

                  self.cachedData = True

      def computeGradientFV3(self, varField):

            runInParallel = False
            ## If we have not actually precomputed necessary datastructures already, let us do it now
            if runInParallel:
                  if not self.cachedData:
                        precomputeGradientFV3Data()
            else:
                  assert(self.cachedData == True)

            # Gradients are 3 component vectors
            nc = 3

            NC = len(self.areaInvD)
            if runInParallel:
                  computepool = multiprocessing.Pool(processes=1)#self.NPROCS)
                  results = computepool.starmap(self.optimizedComputeGradientFV3_Private, 
                                                zip(repeat(varField), self.areaInvD, self.sid, self.fluxIntegral))
                  computepool.close()
                  computepool.join()
                  varGradient = np.reshape(np.array(results, dtype=SF).T, (nc, varField.shape[0]))
            else:
                  varGradient = np.zeros((nc, NC))
                  t_fluxIntegralTotal = np.zeros(nc, dtype=np.float64)
                  gradientLoop(varGradient,t_fluxIntegralTotal, varField, self.areaInvD, self.NS, self.sid, self.fluxIntegral)

            return varGradient

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
            self.areaInvD = []
            self.sid = []
            self.fluxIntegral = []
            self.cachedData = False

      def precomputeGradientFV3Data_Private(self, NC, NP, jj, varCon, varCoords, varStenDex, radius, cellCoords):
            SF = np.float64

            # Gradients are 3 component vectors
            nc = 3

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
            dualEdgeMap = np.zeros((nc,NS))
            boundaryNorm = np.zeros((nc,NS))
            boundaryAngles = np.zeros((NS,1))
            sid = np.zeros((NS,2))
            fluxIntegral = np.zeros((nc,NS), dtype=SF)
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

                  sid[pp, 0] = sid1
                  sid[pp, 1] = sid2

                  fluxIntegral[:,pp] =  vWeight * boundaryNorm[:,pp]

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

      def precomputeGradientFV3Data(self):

            if not self.cachedData:
                  SF = np.float64
                  NC = int(self.varStenDex.shape[0])
                  NP = int(self.varStenDex.shape[1])

                  pool = multiprocessing.Pool(processes=self.NPROCS)
                  results = pool.starmap(self.precomputeGradientFV3Data_Private, zip(repeat(NC), repeat(NP), range(NC), repeat(self.varCon), repeat(self.varCoords), repeat(self.varStenDex), repeat(self.radius), repeat(self.cellCoords)))
                  pool.close()
                  pool.join()

                  ## Store the data in the local object for future use
                  self.areaInvD = np.zeros(NC, dtype=SF)
                  for elem in range(NC):
                        self.areaInvD[elem] = results[elem][0]
                        self.sid.append(np.array(results[elem][1], dtype='i4'))
                        self.fluxIntegral.append(results[elem][2])
                  self.cachedData = True

      def optimizedComputeGradientFV3_Private(self, varField, areaInvD, sid, fluxIntegral):
            SF = np.float64

            # Gradients are 3 component vectors
            nc = 3

            # Loop over the convex hull stencil and get dual edges map
            NS = len(sid[:,0])
            fluxIntegralTotal = np.zeros(nc, dtype=SF)
            for pp in range(NS):

                  # Compute the weighted average of the two cell values AT the shared edge location
                  # VSM: we have precomputed and stored values of sid already
                  varAvg = varField[sid[pp,0]] + varField[sid[pp,1]]

                  # Compute the integral over this edge as a weighted addition with precomputed
                  # fluxIntegral on this edge
                  fluxIntegralTotal = np.add(fluxIntegralTotal, varAvg * fluxIntegral[:,pp])

            # Compute the local gradient at this cell
            varGradient = areaInvD * fluxIntegralTotal

            return varGradient

      def computeGradientFV3(self, varField):

            ## If we have not actually precomputed necessary datastructures already, let us do it now
            if not self.cachedData:
                  self.precomputeGradientFV3Data()

            SF = np.float64

            # Gradients are 3 component vectors
            nc = 3

            computepool = multiprocessing.Pool(processes=self.NPROCS)
            results = computepool.starmap(self.optimizedComputeGradientFV3_Private, 
                                           zip(repeat(varField), self.areaInvD, self.sid, self.fluxIntegral))
            computepool.close()
            computepool.join()
            varGradient = np.reshape(np.array(results, dtype=SF).T, (nc, varField.shape[0]))

            return varGradient

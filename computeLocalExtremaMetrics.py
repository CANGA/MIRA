#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 08:55:56 2018

Computes the local extrema metrics for regridded and reference target data 

@author: jeguerra
"""

import math as mt
import numpy as np
from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

def computeCentroid(NP, cell):
       centroid = np.mat([0.0, 0.0, 0.0])
       for pp in range(NP):
              centroid += cell[:,pp]
              
       centroid *= 1.0 / NP
       
       return centroid

def computeMaximumDistance(NP, pcloud, centroid):
       
       maxDist = 0.0
       cdiff = pcloud - centroid.T
       for pp in range(NP):
              dist = mt.sqrt(cdiff[0,pp]**2 + cdiff[1,pp]**2 + cdiff[2,pp]**2)
              
              if dist > maxDist:
                     maxDist = dist
                     
       return maxDist

def computeCoordPatchIndexArray(NC, pcloud, centroid, radius):
       
       pdex = []
       cdiff = pcloud - centroid.T
       for ii in range(NC):
              dist = mt.sqrt(cdiff[0,ii]**2 + cdiff[1,ii]**2 + cdiff[2,ii]**2)
              
              if dist < radius:
                     pdex.append(ii)
       
       return pdex

def computeLocalPatchExtrema(jj, varConS, varCoordS, varS, varConT, varCoordT):
       
       # Index the target cell from the input coordinates
       cellT = varConT[jj,:]
       NP = len(cellT)
       cdexT = cellT - 1;
       cell = varCoordT[:,cdexT]
       
       # compute the centroid of the target cell
       centroid = computeCentroid(NP, cell)
       
       # Compute the maximum distance from the centroid to the corners of the target cell
       mdistT = computeMaximumDistance(NP, cell, centroid)
                     
       # Search the source coordinates for any that are within the radius mdistT
       NS = varCoordS.shape[1]
       odexS = computeCoordPatchIndexArray(NS, varCoordS, centroid, mdistT)
                     
       # Find the source cells (patch) that have any of the coordinates found above
       NS = varConS.shape[0]
       cdexS = []
       for oo in range(len(odexS)):
              for ii in range(NS):
                     # Use only the first connectivity column
                     if varConS[ii,0] == odexS[oo]:
                            cdexS.append(ii)
                            
       # Compute local max and min from source data in the patch defined above
       pmax = np.amax(varS[cdexS])
       pmin = np.amin(varS[cdexS])
              
       return pmin, pmax

def computeLocalExtremaMetrics(areaT, varSS, varS2T, varST, varConS, varCoordS, varConT, varCoordT):
       
       NT = len(varST)
       minDiff = np.zeros((NT,1))
       maxDiff = np.zeros((NT,1))
       # Compute the localized difference arrays (eqs. 10 and 11)
       #for jj in range(NT):
       for jj in range(1):
              
              # Compute the patch extrema (expensive calculation)
              lPmin, lPmax = computeLocalPatchExtrema(jj, varConS, varCoordS, varSS, varConT, varCoordT)
              
              # Compute the min and max difference arrays
              minDiff[jj] = np.minimum(varST[jj] - lPmin, 0.0)
              maxDiff[jj] = np.maximum(varST[jj] - lPmax, 0.0)
              
       # Compute standard norms on local extrema differences
       NT = len(varST)
       varST2 = np.power(varST, 2)
       
       # Compute normalization integrals
       L1Den = computeGlobalWeightedIntegral(NT, varST, areaT)
       L2Den = computeGlobalWeightedIntegral(NT, varST2, areaT)
       LinfDen = np.amax(varST)
       
       # Compute numerators for minima
       varDiff = minDiff
       varDiff2 = np.power(minDiff, 2)
       L1Num = computeGlobalWeightedIntegral(NT, varDiff, areaT)
       L2Num = computeGlobalWeightedIntegral(NT, varDiff2, areaT)
       LinfNum = np.amax(abs(varDiff))
       
       Lmin_1 = float(L1Num / L1Den)
       Lmin_2 = mt.sqrt(L2Num) / mt.sqrt(L2Den)
       Lmin_inf = LinfNum / LinfDen
       
       # Compute numerators for maxima
       varDiff = maxDiff
       varDiff2 = np.power(maxDiff, 2)
       L1Num = computeGlobalWeightedIntegral(NT, varDiff, areaT)
       L2Num = computeGlobalWeightedIntegral(NT, varDiff2, areaT)
       LinfNum = np.amax(abs(varDiff))
       
       Lmax_1 = float(L1Num / L1Den)
       Lmax_2 = mt.sqrt(L2Num) / mt.sqrt(L2Den)
       Lmax_inf = LinfNum / LinfDen
       
       return Lmin_1, Lmin_2, Lmin_inf, Lmax_1, Lmax_2, Lmax_inf 
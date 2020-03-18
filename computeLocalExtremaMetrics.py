#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 08:55:56 2018

Computes the local extrema metrics for regridded and reference target data 

@author: jeguerra
"""

import math as mt
import numpy as np
#from scipy.spatial import cKDTree
from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

COINCIDENT_TOLERANCE = 1.0E-14
kdleafs = 100

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

''' OLD METHOD INVOLVING A PATCH SEARCH OF THE SOURCE MESH
def computeLocalPatchExtrema(jj, varConS, coordTree, varS, varConT, varCoordT):
       
       # Index the target cell from the input coordinates
       cellT = varConT[jj,:]
       NP = len(cellT)
       cdexT = cellT - 1;
       cell = varCoordT[:,cdexT]
       
       # compute the centroid of the target cell
       centroid = computeCentroid(NP, cell)
       
       # Compute the maximum distance from the centroid to the corners of the target cell
       mdistT = computeMaximumDistance(NP, cell, centroid)
                     
       # Compute the patch of nearest nodes to centroid in the source mesh
       ndex = coordTree.query_ball_point(centroid, mdistT, p=2, eps=0)
       ndex = ndex[0]
       # Compute the patch of nearest cells based on previous result
       NN = len(ndex)
       nn = np.add(ndex,1)
       
       cdex = []
       # Loop over any nodes found to be overlapping jj cell
       for ii in range(NN):
              # Loop over each column of the nodal connectivity
              for cc in range(NP):
                     idex = np.where(varConS[:,cc] == nn[ii])
                     idex = np.ravel(idex[0])
                     
                     if len(idex) == 0:
                            continue
                     else:
                            cdex.append(idex.tolist())
       
       # Fancy Python bit to merge the patch indices for source cells
       cdex = [item for sublist in cdex for item in sublist]
       
       if len(cdex) >= 1:
              pmin = np.amin(varS[cdex])
              pmax = np.amax(varS[cdex])
       else:
              # In case no patch was found
              pmin = 0.0
              pmax = 0.0
              
       return pmin, pmax
'''

def computeLocalPatchExtrema(jj, varConStenDexT, varConT, varST, SpectralElement):
       
       # Fetch the stencil of neighboring elements/cells
       sdex = varConStenDexT[jj,:] - 1
       sdex = sdex.astype(int)

       if SpectralElement:
              # Fetch the gridID for all nodes in the stencil of elements
              ndex = varConT[sdex,:] - 1
              ndex = ndex.astype(int)
              # Fetch the nodal values 
              varPatch = varST[ndex]
       else:
              # Fetch the cell values
              varPatch = varST[sdex]

       pmin = np.amin(varPatch)
       pmax = np.amax(varPatch)
       
       return pmin, pmax

def computeLocalExtremaMetrics(varConStenDex, varCon, varCoord, varS2T, varST, areaT, jacobiansT, SpectralElement):
       
       NT = varCon.shape[0]
       varST2 = np.power(varST, 2)
       minDiff = np.zeros((NT,1))
       maxDiff = np.zeros((NT,1))
       
       # Compute a KDtree for the source coordinates
       #coordTreeS = cKDTree(varCoordS.T, leafsize=kdleafs)
       
       # Compute the localized difference arrays (eqs. 10 and 11)
       for jj in range(NT):
              
              # Compute the patch extrema using the KDtree set up above (OLD WAY USING SOURCE DATA)
              #lPmin, lPmax = computeLocalPatchExtrema(jj, varConS, coordTreeS, varSS, varConT, varCoordT)
              
              # Compute the patch extrema using the sampled target data and adjacency stencil
              lPmin, lPmax = computeLocalPatchExtrema(jj, varConStenDex, varCon, varST, SpectralElement)

              # Compute the min and max difference arrays
              minDiff[jj] = np.abs(varS2T[jj] - lPmin)
              maxDiff[jj] = np.abs(lPmax - varS2T[jj])
       
       # Compute normalization integrals
       L1Den = computeGlobalWeightedIntegral(NT, varCon, np.abs(varST), areaT, jacobiansT, SpectralElement)
       L2Den = computeGlobalWeightedIntegral(NT, varCon, np.abs(varST2), areaT, jacobiansT, SpectralElement)
       LinfDen = np.amax(varST) - np.amin(varST)
       
       # Compute numerators for minima
       varDiff = minDiff
       varDiff2 = np.power(minDiff, 2)
       L1Num = computeGlobalWeightedIntegral(NT, varCon, varDiff, areaT, jacobiansT, SpectralElement)
       L2Num = computeGlobalWeightedIntegral(NT, varCon, varDiff2, areaT, jacobiansT, SpectralElement)
       LinfNum = np.amax(varDiff)

       LareaDen = computeGlobalWeightedIntegral(NT, varCon, np.ones(varST.shape), areaT, jacobiansT, SpectralElement)
       
       if SpectralElement:
              Lmin_1 = np.asscalar(L1Num / L1Den)
              Lmin_2 = mt.sqrt(L2Num / L2Den)
              Lmin_inf = LinfNum / LinfDen
       else:
              Lmin_1 = np.asscalar(L1Num / NT / LinfDen / LareaDen)
              Lmin_2 = mt.sqrt(L2Num / NT) / LinfDen / LareaDen
              Lmin_inf = LinfNum / LinfDen

       # Compute numerators for maxima
       varDiff = maxDiff
       varDiff2 = np.power(maxDiff, 2)
       L1Num = computeGlobalWeightedIntegral(NT, varCon, varDiff, areaT, jacobiansT, SpectralElement)
       L2Num = computeGlobalWeightedIntegral(NT, varCon, varDiff2, areaT, jacobiansT, SpectralElement)
       LinfNum = np.amax(varDiff)
       
       if SpectralElement:
              Lmax_1 = np.asscalar(L1Num / L1Den)
              Lmax_2 = mt.sqrt(L2Num / L2Den)
              Lmax_inf = LinfNum / LinfDen
       else:
              Lmax_1 = np.asscalar(L1Num / NT / LinfDen / LareaDen)
              Lmax_2 = mt.sqrt(L2Num / NT) / LinfDen / LareaDen
              Lmax_inf = LinfNum / LinfDen
       
       return Lmin_1, Lmin_2, Lmin_inf, Lmax_1, Lmax_2, Lmax_inf
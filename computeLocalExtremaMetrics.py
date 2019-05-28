#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 08:55:56 2018

Computes the local extrema metrics for regridded and reference target data 

@author: jeguerra
"""

import math as mt
import numpy as np
from scipy.spatial import cKDTree
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
              varPatch = varST[ndex.flatten]
       else:
              # Fetch the cell values
              varPatch = varST[sdex]
              
       pmin = np.amin(varPatch)
       pmax = np.amin(varPatch)
       
       return pmin, pmax

def computeLocalExtremaMetrics(varConStenDexT, varConT, varCoordT, varS2T, varST, SpectralElement):
       
       NT = len(varST)
       minDiff = np.zeros((NT,1))
       maxDiff = np.zeros((NT,1))
       
       # Compute a KDtree for the source coordinates
       #coordTreeS = cKDTree(varCoordS.T, leafsize=kdleafs)
       
       # Compute the localized difference arrays (eqs. 10 and 11)
       for jj in range(NT):
              
              # Compute the patch extrema using the KDtree set up above (OLD WAY USING SOURCE DATA)
              #lPmin, lPmax = computeLocalPatchExtrema(jj, varConS, coordTreeS, varSS, varConT, varCoordT)
              
              # Compute the patch extrema using the sampled target data and adjacency stencil
              lPmin, lPmax = computeLocalPatchExtrema(jj, varConStenDexT, varConT, varST, SpectralElement)
              
              # Compute the min and max difference arrays
              minDiff[jj] = np.minimum(varS2T[jj] - lPmin, 0.0)
              maxDiff[jj] = np.maximum(lPmax - varS2T[jj], 0.0)
              
       # Compute standard norms on local extrema differences
       NT = len(varST)
       
       # Compute normalization integrals
       LinfDen = np.amax(varST)
       
       # Compute numerators for minima
       varDiff = minDiff
       LinfNum = np.amax(abs(varDiff))
       
       Lmin_inf = LinfNum / LinfDen
       
       # Compute numerators for maxima
       varDiff = maxDiff
       LinfNum = np.amax(abs(varDiff))
       
       Lmax_inf = LinfNum / LinfDen
       
       return Lmin_inf, Lmax_inf 
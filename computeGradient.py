#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:26:47 2018

Compute a low order FV gradient on the manifold based on the local intersection
of the planes defined by two connected cells sharing an edge. plane1: defined
by the centroids of the cells and plane2: defined by the nodes on the shared edge

By construction, all the vectors involved start at the origin... convenient.
Parameterized flux integral assumes constant radius and field value.

@author: jeguerra
"""

import numpy as np
import math as mt
from computeEdgesArray import computeEdgesArray

def computeCentroid(NP, cell):
       centroid = np.mat([0.0, 0.0, 0.0])
       for pp in range(NP):
              centroid += cell[:,pp]
              
       centroid *= 1.0 / NP
       
       return centroid

def computeGradient(varList, varCoords, varStenDex, areas):
       
       varGradient = [np.zeros(np.size(varList[0])), \
                      np.zeros(np.size(varList[1]))]
       
       cellCoords = np.zeros((3,areas.shape[0]))
       
       NV = len(varList)
       NC = varStenDex.shape[0]
       NP = int(varStenDex.shape[1] / 2)
       # Loop over the cells
       for jj in range(NC):
              pdex = np.array(range(NP), dtype = int)
              # Compute the center cell centroid
              cdex = (varStenDex[jj, pdex]) - 1
              cdex = cdex.astype(int)
              cell = varCoords[:,cdex]
              centroidC = computeCentroid(NP, cell)
              radiusC = np.linalg.norm(centroidC)
              
              # Set the centroid coordinates to the workspace
              cellCoords[:,jj] = centroidC
              
              # Get the local node pair map for these edges (indices)
              edgeDex = computeEdgesArray(NP, (cdex + 1))
              
              # Check for local degeneracy in stencil and fix connectivity
              for pp in range(NP):
                     # Look for 0 in the adjacency stencil
                     if varStenDex[jj,NP + pp] == 0:
                            pdex = np.delete(pdex, pp)
                     else:
                            continue
              
              # Loop over the stencil of the gradient
              fluxIntegral = np.zeros(NV)
              for pp in pdex:
                     # Get the ID and centroid of the connected cell
                     sid = varStenDex[jj, NP + pp] - 1
                     sid = sid.astype(int)
                     cdex = (varStenDex[sid, range(NP)]) - 1
                     cdex = cdex.astype(int)
                     cell = varCoords[:,cdex]
                     centroidS = computeCentroid(NP, cell)
                     radiusS = np.linalg.norm(centroidS)
                     
                     # Get coordinates of the shared edge nodes
                     nid1 = edgeDex[pp,0] - 1
                     nid2 = edgeDex[pp,1] - 1
                     coord1 = varCoords[:, nid1]
                     coord2 = varCoords[:, nid2]
                     
                     # Compute normals of the intersecting planes
                     cellCross = np.cross(centroidC, centroidS)
                     nodeCross = np.cross(coord1, coord2)
                     """
                     print(coord1)
                     print(' ')
                     print(coord2)
                     print(' ')
                     print(nodeCross)
                     print('--------')
                     print(centroidC)
                     print(' ')
                     print(centroidS)
                     print(' ')
                     print(cellCross)
                     print('********')
                     """
                     # Normalize the... normal vectors =)
                     unCell = 1.0 / np.linalg.norm(cellCross) * cellCross
                     unNode = 1.0 / np.linalg.norm(nodeCross) * nodeCross
                     
                     # Compute unit vector that intersects the shared edge
                     unEdge = np.cross(unCell, unNode)
                     
                     # After all that geometric juggling...
                     
                     # Get the angle spanned by going from coord1 to coord2 and average local radius
                     RE = 0.5 * (np.linalg.norm(coord1) + np.linalg.norm(coord2))
                     unCord1 = 1.0 / np.linalg.norm(coord1) * coord1
                     unCord2 = 1.0 / np.linalg.norm(coord2) * coord2
                     Alpha = mt.acos(np.dot(unCord1, unCord2))
                     Alpha = abs(Alpha)
                     
                     # Get the angle spanned by going from target cell centroid to unEdge
                     unCent = 1.0 / np.linalg.norm(centroidC) * centroidC
                     vdot = np.dot(np.ravel(unCent), np.ravel(unEdge))
                     beta = mt.acos(vdot)
                     beta = abs(beta)
                     
                     # Get the total angle separating the two cell centroids
                     #print(jj, radiusC, radiusS)
                     vdot = np.dot(np.ravel(centroidC), np.ravel(centroidS))
                     Beta = mt.acos(vdot / radiusC / radiusS)
                     Beta = abs(Beta)
                     
                     # Check alpha for zero or < zero, make positive definite
                     if abs(Alpha) < 1.0E-6:
                            print('Sliver cell detected! Check your mesh/code at cell: ', (jj + 1))
                            
                     for vv in range(NV):
                            # Compute the weighted average of the two cell values AT the shared edge location
                            vWeight = beta / Beta
                            varAvg = varList[vv][jj] + \
                                     vWeight * (varList[vv][sid] - varList[vv][jj])
                            
                            # Compute the integral over this edge
                            fluxIntegral[vv] += RE * Alpha * varAvg
                            
              
              # Compute the local gradient at this cell
              for vv in range(NV):
                     varGradient[vv][jj] = 1.0 / areas[jj] * fluxIntegral[vv]
                     #print('Variable: ', varList[vv][jj])
                     #print('Gradient: ', varGradient[vv][jj])
              
       return varGradient, cellCoords

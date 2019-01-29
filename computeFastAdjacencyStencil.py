#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 21 09:16:45 2018

Computes the adjacency array (augmented to the connectivity). Assumes that all
cells have the same orientation (clockwise or aniclockwise). For a given cell,
make an (NP, 2) array of edges in reversed order to match any other cells that
also have that same edge, again assuming all cells have the same orientation.
Supports: quadrilaterals and triangles.

NOTE: edgeMap: [cell 2, node1, node2, cell 1] outer give the two cells that
belong to the edge specified in nodes. left to right connectivity for cell 1
and right to left connectivity for cell 2 

@author: jeguerra
"""

import numpy as np
from scipy.spatial import cKDTree
from computeEdgesArray import computeEdgesArray

def computeFastAdjacencyStencil(varCon):
       
       COINCIDENT_TOLERANCE = 1.0E-14
       NC = varCon.shape[0]
       NP = varCon.shape[1]
       
       varConStenDex = np.zeros((NC, NP + NP))
       
       # Make an array of edges based on grid pairs from connectivity and a cell id
       # This has coincident pairs of edges for each cell processed
       edgeMap = np.zeros((NP, 3))
       for cc in range(NC):
              
              # Copy over the connectivity of the current cell
              varConStenDex[cc, range(NP)] = varCon[cc,:]
              
              # Make the cell id column for dim 3 of edgeMap
              cid = (cc + 1) * np.ones((NP, 1))
              
              # Get the local node pair map for these edges
              edges = computeEdgesArray(NP, varCon[cc,:])
              
              # Append the cell map to the end of the node map
              edges = np.append(edges, cid, axis=1)
              
              if cc == 0:
                     edgeNodeMap = edges
              else:
                     edgeNodeMap = np.append(edgeNodeMap, edges, axis=0)
                     
       # Sort edgeMap first 2 column in ascending order to enable coincident check
       sortedEdges = np.sort(edgeNodeMap[:,[0, 1]], axis=1)
                     
       # Build a KDtree around the edge map of cell - node pair
       edgeTree = cKDTree(sortedEdges, leafsize=100)
       
       # Compute the edges to cells map [c1 n1 n2 c2]
       NE = np.size(edgeNodeMap, axis=0)
       for ii in range(NE):
              # Compute the list of nodes that are coincident
              thisEdge = sortedEdges[ii,:]
              
              # ndex will have 2 cellID's for non-degenerate edges
              ndex = edgeTree.query_ball_point(thisEdge, COINCIDENT_TOLERANCE, p=2, eps=0)
              
              # coinDex stores indices of coincident edges (and degenerate points)
              # keepDex stores indices for unique edges
              if ii == 0:
                     coinDex = [ndex]
                     keepDex = [int(min(ndex))]
                     edgeCellMap = np.array([[edgeNodeMap[ndex[0],2], \
                                              thisEdge[0], thisEdge[1], \
                                              edgeNodeMap[ndex[1],2]]])
              else:
                     coinDex.append(ndex)
                     keepDex.append(int(min(ndex)))
                     thisCellMap = np.array([edgeNodeMap[ndex[0],2], \
                                             thisEdge[0], thisEdge[1], \
                                             edgeNodeMap[ndex[1],2]])
                     edgeCellMap = np.append(edgeCellMap, [thisCellMap], axis=0)
                     
       # Trim coincidents from edgeCellMap with keepDex array
       keepDex = np.unique(keepDex)
       cleanEdgeCellMap = edgeCellMap[keepDex,:]
                     
       # Compute a KDtree from the smaller cleanEdgeCellMap
       edgeTree = cKDTree(cleanEdgeCellMap, leafsize=100)
       
       # Loop over the node connectivity and construct the adjacency stencil
                            
       return edgeNodeMap, edgeCellMap, cleanEdgeCellMap, coinDex, varConStenDex

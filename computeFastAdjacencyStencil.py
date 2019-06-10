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

COINCIDENT_TOLERANCE = 1.0E-14
kdleafs = 64

def computeFastAdjacencyStencil(varCon):
       
       NC = varCon.shape[0]
       NP = varCon.shape[1]
       
       varConStenDex = np.zeros((NC, NP + NP))
       
       # Make an array of edges based on grid pairs from connectivity and a cell id
       # This has coincident pairs of edges for each cell processed
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
                     
       edgeTree = cKDTree(edgeNodeMap[:,[0, 1]], leafsize=kdleafs)
       # Loop over the node connectivity and construct the adjacency stencil
       for ii in range(NC):
              # Get the local node pair map for these edges
              edges = computeEdgesArray(NP, varCon[ii,:])
              
              # Loop over the surrounding edges to this cell
              for jj in range(NP):
                     # Check for degenerate edge leaves a 0 in the stencil
                     if edges[jj,0] == edges[jj,1]:
                            continue
                     
                     # Fetch the current edge in both local directions
                     thisEdge = edges[jj,::-1]
                     
                     # Find the matching edge (should only give one result)
                     cdex = edgeTree.query_ball_point(thisEdge, COINCIDENT_TOLERANCE, p=2, eps=0)
                     
                     # Check for no edge found (indicates a hole)
                     if not cdex:
                            continue
                     
                     # Get the connected cell and set the stencil
                     varConStenDex[ii,NP+jj] = edgeNodeMap[cdex,2]
                            
       return edgeNodeMap, edgeTree, varConStenDex

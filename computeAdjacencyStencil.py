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

def computeEdgesArray(NP, connect, oppose):
       
       if NP == 4:
              # Set the direction of edges
              if oppose == True:
                     edex = [1, 0, 2, 1, 3, 2, 0, 3]
              elif oppose == False:
                     edex = [0, 1, 1, 2, 2, 3, 3, 0]
              else:
                     edex = [0, 1, 1, 2, 2, 3, 3, 0]
              
              edges = np.array([[connect.flat[edex[0]], connect.flat[edex[1]]], \
                                [connect.flat[edex[2]], connect.flat[edex[3]]], \
                                [connect.flat[edex[4]], connect.flat[edex[5]]], \
                                [connect.flat[edex[6]], connect.flat[edex[7]]]])
       elif NP == 3:
              # Set the direction of edges
              if oppose == True:
                     edex = [1, 0, 2, 1, 0, 2]
              elif oppose == False:
                     edex = [0, 1, 1, 2, 2, 0]
              else:
                     edex = [0, 1, 1, 2, 2, 0]
                     
              edges = np.array([[connect.flat[edex[0]], connect.flat[edex[1]]], \
                                [connect.flat[edex[2]], connect.flat[edex[3]]], \
                                [connect.flat[edex[4]], connect.flat[edex[5]]]])
       else:
              print('ONLY TRIANGLES OR QUADRILATERALS SUPPORTED!')
       
       return edges

def computeAdjacencyStencil(varCon):
       
       NC = varCon.shape[0]
       NP = varCon.shape[1]
       
       varConStenDex = np.zeros((NC, NP + NP))
       
       # Make an array of edges based on grid pairs from connectivity and a cell id
       # This has coincident pairs of edges for each cell processed
       edgeMap = np.zeros((NP, 3))
       secondCell = []
       for cc in range(NC):
              
              # Copy over the connectivity of the current cell
              varConStenDex[cc, range(NP)] = varCon[cc,:]
              
              # Make the cell id column for dim 3 of edgeMap
              cid = (cc + 1) * np.ones((NP, 1))
              
              # Get the local node pair map for these edges
              edges = computeEdgesArray(NP, varCon[cc,:], False)
              
              # Append the cell map to the end of the node map
              edges = np.append(edges, cid, axis=1)
              
              if cc == 0:
                     edgeMap = edges
              else:
                     edgeMap = np.append(edgeMap, edges, axis=0)
                     
       # Loop over the full edge map with coincident edges
       # Make the edges to cells map AND remove coincident edges
       NE = edgeMap.shape[0]
       NI = NE
       dd = 0
       while dd < NE / 2:
              # Fetch the current edge
              refEdge = edgeMap[dd,[0,1]]
              
              # Loop over modified edgeMap
              for ee in range(dd+1, NI):
                     
                     # Fetch a comparison edge
                     cmpEdge = edgeMap[ee,[0,1]]
                     # Make the comparison of the node map for coincident edges
                     n1Diff = abs(cmpEdge[0] - refEdge[1])
                     n2Diff = abs(cmpEdge[1] - refEdge[0])
                     
                     if (n1Diff == 0) and (n2Diff == 0):
                            # Get the other cell that belongs to refEdge
                            secondCell = np.append(secondCell, edgeMap[ee,2])
                            # Remove the coincident edge
                            edgeMap = np.delete(edgeMap, ee, axis=0)
                            NI = edgeMap.shape[0]
                            dd += 1
                            break
       
       # Append the second cell to the edgeMap
       edgeMap = np.insert(edgeMap, 0, secondCell.T, axis=1)
       
       # Loop over the right most column of the edgeMap and set stencil partially
       NM = edgeMap.shape[0]
       mdex = np.zeros(NC, dtype=int)
       kk = 0
       cdex = 0
       for mm in range(NM):
              
              if mm > 0:
                     # Check that we are on the same cell else reset kk
                     if int(edgeMap[mm,3]) - 1 != cdex:
                            kk = 0
              
              # Get index of the current cell (right column of edgeMap)
              cdex = int(edgeMap[mm,3]) - 1
              # Get index of the connected cell (left column of edgeMap)
              sdex = int(edgeMap[mm,0])
              
              varConStenDex[cdex, NP + kk] = sdex
              kk += 1
              mdex[cdex] = kk
              
       # Fetch the left column of edgeMap and argsort() to get sorting indices
       missingCells = edgeMap[:, 0]
       sortDex = np.argsort(missingCells)
       
       # Apply sorting to edgeMap -> sortedEdgeMap
       sortedEdgeMap = edgeMap[sortDex,:]
       
       # Loop over the sorted edgeMap using the left column to set the remainder stencil
       kk = 0
       cdex = 0
       for mm in range(NM):
              
              if mm > 0:
                     # Check that we are on the same cell else reset kk
                     if int(edgeMap[mm,0]) - 1 != cdex:
                            kk = 0
              
              # Get index of the current cell (right column of edgeMap)
              cdex = int(sortedEdgeMap[mm,0]) - 1
              # Get index of the connected cell (left column of edgeMap)
              sdex = int(sortedEdgeMap[mm,3])
              
              edex = (NP + mdex[cdex] + kk)
              varConStenDex[cdex, edex] = sdex
              kk += 1
              """
              if (mdex[cdex] > (NP - 1)) and (mdex[cdex] < 2 * NP):
                     varConStenDex[cdex, mdex[cdex] + kk] = sdex
              else:
                     continue
              """
                            
       return mdex, edgeMap, sortedEdgeMap, varConStenDex

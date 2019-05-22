#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 11:05:32 2019

Computes new global GLL coordinates and connectivity using edge map from adjacency

@author: jeguerra
"""

import numpy as np
import math as mt
from scipy import optimize

COINCIDENT_TOLERANCE = 1.0E-14

def computeEdgeParametersGLL(coord1, coord2):
       # Define arc segment from thisEdge[jj,0] to thisEdge[jj,1]
       RE1 = np.linalg.norm(coord1)
       RE2 = np.linalg.norm(coord2)
       RE = 0.5 * (RE1 + RE2)
       # Compute unit position vectors for the edge end points
       unCoord1 = 1.0 / RE1 * coord1
       unCoord2 = 1.0 / RE2 * coord2
       edgeAngle = mt.acos(np.dot(unCoord1, unCoord2))
       
       # Compute the plane where the edge lies
       edgePlane = np.cross(unCoord1, unCoord2)
       
       # Compute the half way coordinate (for initial guess)
       halfCoord = 0.5 * (coord1 + coord2)
       halfCoord *= RE / np.linalg.norm(halfCoord);
       
       return halfCoord, edgePlane, edgeAngle, RE

def findCoordAlongEdge(x, edgePlane, x1Coord, locAngle, RE):
       
       # Equation 1 = the angle between x and x1Coord is locAngle
       # Equation 2 = the new grid must be on the same plane as the edge
       # Equation 3 = the point is on a sphere of radius RE
       
       return [x1Coord[0] * x[0] + x1Coord[1] * x[1] + x1Coord[2] * x[2] - mt.cos(locAngle), \
               (edgePlane[1] * x1Coord[2] - edgePlane[2] * x1Coord[1]) * x[0] + \
               (edgePlane[2] * x1Coord[0] - edgePlane[0] * x1Coord[2]) * x[1] + \
               (edgePlane[0] * x1Coord[1] - edgePlane[1] * x1Coord[0]) * x[2], \
               x[0]**2 + x[1]**2 + x[2]**2 - RE]

def computeCoordConnGLL(NEL, NGED, NGEL, varCoord, varCon, edgeNodeMap, edgeNodeKDTree, seOrder):
       
       #NEL = total number of elements
       #NGED = number of grids per edge (new mesh)
       #NGEL = number of grids per element (new mesh)
       #NNG = total number of new grids (besides corners)
       
       # Get the number of grids in original mesh
       NG = varCoord.shape[1]
       # Get the number of grids/edges per element in original mesh
       NEEL = varCon.shape[1]
       # Get the new number of nodes per element on edges only
       NGEO = NEEL * (NGED - 1)
       
       #""" 2nd order method
       if seOrder == 2:
              GN = [-1.0, \
                     0.0, \
                    +1.0]
              
       #""" 4th order method
       if seOrder == 4:
              GN = [-1.0, \
                    -mt.sqrt(1.0 / 5.0), \
                    +mt.sqrt(1.0 / 5.0), \
                    +1.0]
              
       # Change GN to [0.0 1.0]
       GN = 0.5 * np.add(GN, 1.0)
       
       # Initialize new global GLL connectivity
       varConGLL = np.zeros((NEL,NGEL))
       # Initialize new global GLL complement of grids
       varCoordGLL = varCoord
       # Initialize new global GLL edge-node map (last column is the cell id)
       NED = edgeNodeMap.shape[0]
       edgeNodeMapGLL = np.zeros((NED,NGED))
              
       # Loop over the edges
       gg = 0
       for ii in range(NED):
              # Get the local node pair map for these edges
              thisEdge = edgeNodeMap[ii,0:2]
              
              # Get the opposing edge from neighbor element
              thatEdge = thisEdge[::-1]
              
              # Set the end points to the new GLL connectivity
              edgeNodeMapGLL[ii,0] = int(thisEdge[0])
              edgeNodeMapGLL[ii,NGED-1] = int(thisEdge[1])
              
              # Check for degenerate edge leaves a 0 in the connectivity
              if thisEdge[0] == thisEdge[1]:
                     continue
              
              # Check that edgeNodeMapGLL hasn't been already set for this edge
              if sum(edgeNodeMapGLL[ii,1:NGED-1]) != 0.0:
                     continue
              
              # Fetch the end point coordinates
              coord1 = varCoord[:,int(thisEdge[0])-1]
              coord2 = varCoord[:,int(thisEdge[1])-1]
              
              halfCoord, edgePlane, edgeAngle, RE = computeEdgeParametersGLL(coord1, coord2)
              
              # Loop over the new grids on the edge
              for jj in range(1,NGED-1):
                     # Compute angular location of new grid on the edge
                     newGridAngle = GN[jj] * edgeAngle
                     # Solve for the new coordinate
                     sol = optimize.root(findCoordAlongEdge, halfCoord, (edgePlane, coord1, newGridAngle, RE))
                     # Store the new grid
                     newGrid = np.zeros((3,1))
                     newGrid[:,0] = sol.x
                     varCoordGLL = np.append(varCoordGLL, newGrid, axis=1)
                     # Make a new grid ID
                     newGridID = NG + gg + 1
                     # Update the new connectivity
                     edgeNodeMapGLL[ii,jj] = int(newGridID)
                     gg += 1
                     
              # Set connectivity for the opposing edge in edgeNodeMap (thatEdge)
              edex = edgeNodeKDTree.query_ball_point(thatEdge, COINCIDENT_TOLERANCE, p=2, eps=0)
              edgeNodeMapGLL[edex,1:NGED-1] = edgeNodeMapGLL[ii,NGED-2:0:-1]
              
       print(varCoord.shape, varCoordGLL.shape)       
       # Loop over the elements and reconstruct new connectivity for edges only
       gg = 0
       edex = range(NEEL)
       for ii in range(NEL):
              # Fectch every NEEL edges per element
              if ii > 0:
                     edex = np.add(edex, NEEL)
              thisElement = edgeNodeMapGLL[edex,:]
              
              # Loop over the perimeter edges setting new interior nodes to the connectivity
              cdex =range(1,NEEL)
              for jj in range(NEEL):
                     if jj == 0:
                            varConGLL[ii,0:NEEL] = thisElement[0,:]
                     elif jj == NEEL - 1:
                            varConGLL[ii,(NGEO-2):NGEO] = thisElement[jj,1:NEEL-1]
                     else:
                            cdex = np.add(cdex, (NEEL - 1))
                            varConGLL[ii,cdex] = thisElement[jj,1:NEEL]
                            
              # Now set the element interior GLL nodes using edges 2 and 4
              # Order follows connectivity of the parent element
              edge24 = np.array([[thisElement[3,2], thisElement[1,1]], \
                        [thisElement[1,2], thisElement[3,1]]])
              
              # Fetch the end point coordinates (3, 2, 2) array
              coords = np.array([[varCoordGLL[:,int(edge24[0,0])-1], 
                                  varCoordGLL[:,int(edge24[0,1])-1]], \
                                 [varCoordGLL[:,int(edge24[1,0])-1], 
                                  varCoordGLL[:,int(edge24[1,1])-1]]])
                     
              # Loop over the rows (2nd dim) of coords and make new global nodes
              hh = 0
              for jj in range(coords.shape[1]):
                     coord1 = coords[jj,0,:]
                     coord2 = coords[jj,1,:]
                     halfCoord, edgePlane, edgeAngle, RE = computeEdgeParametersGLL(coord1, coord2)
                                          
                     # Loop over the new grids on the edge
                     for kk in range(1,NGED-1):
                            # Compute angular location of new grid on the edge
                            newGridAngle = GN[jj] * edgeAngle
                            # Solve for the new coordinate
                            sol = optimize.root(findCoordAlongEdge, halfCoord, (edgePlane, coord1, newGridAngle, RE))
                            # Store the new grid
                            newGrid = np.zeros((3,1))
                            newGrid[:,0] = sol.x
                            varCoordGLL = np.append(varCoordGLL, newGrid, axis=1)
                            # Make a new grid ID (continuing from above)
                            newGridID += 1
                            #print(newGridID)
                            #print(newGrid)
                            # Put the new grid ID at the end of the connectivity row
                            varConGLL[ii,NGEO+hh] = newGridID
                            hh += 1
                            gg += 1
                            
                            
       return edgeNodeMapGLL, varCoordGLL, varConGLL

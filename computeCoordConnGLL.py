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

def findCoordAlongEdge(x, edgePlane, x1Coord, locAngle, RE):
       
       # Equation 1 = the angle between x and x1Coord is locAngle
       # Equation 2 = the new grid must be on the same plane as the edge
       # Equation 3 = the point is on a sphere of radius RE
       
       return [x1Coord[0] * x[0] + x1Coord[1] * x[1] + x1Coord[2] * x[2] - mt.cos(locAngle), \
               (edgePlane[1] * x1Coord[2] - edgePlane[2] * x1Coord[1]) * x[0] + \
               (edgePlane[2] * x1Coord[0] - edgePlane[0] * x1Coord[2]) * x[1] + \
               (edgePlane[0] * x1Coord[1] - edgePlane[1] * x1Coord[0]) * x[2], \
               x[0]**2 + x[1]**2 + x[2]**2 - RE]

def computeCoordConnGLL(NEL, NGED, NGEL, NNG, varCoords, edgeNodeMap, edgeNodeKDTree, seOrder):
       
       #NEL = total number of elements
       #NGED = number of grids per edge
       #NGEL = number of grids per element
       #NNG = total number of new grids (besides corners)
       
       NG = varCoords.shape[1]
       
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
       varCoordGLL = np.zeros((3,NNG))
       # Initialize new global GLL edge-node map (last column is the cell id)
       NED = edgeNodeMap.shape[0]
       edgeNodeMapGLL = np.zeros((NED,NGED + 1))
              
       # Loop over the edges
       gg = 0
       #for ii in range(NED):
       for ii in range(10):
              # Get the local node pair map for these edges
              thisEdge = edgeNodeMap[ii,0:2]
              
              # Check for degenerate edge leaves a 0 in the stencil
              if thisEdge[0] == thisEdge[1]:
                     continue
              
              # Fetch the end point coordinates
              coord1 = varCoords[:,int(thisEdge[0])-1]
              coord2 = varCoords[:,int(thisEdge[1])-1]
              
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
              
              # Set the end points to the new GLL connectivity
              edgeNodeMapGLL[ii,0] = int(thisEdge[0])
              edgeNodeMapGLL[ii,NGED-1] = int(thisEdge[1])
              
              # Loop over the new grids on the edge
              for jj in range(1,NGED-1):
                     # Compute angular location of new grid on the edge
                     newGridAngle = GN[jj] * edgeAngle
                     # Solve for the new coordinate
                     sol = optimize.root(findCoordAlongEdge, halfCoord, (edgePlane, coord1, newGridAngle, RE))
                     # Store the new grid
                     varCoordGLL[:,gg] = sol.x
                     # Make a new grid ID
                     newGridID = NG + gg + 1
                     # Update the new connectivity
                     edgeNodeMapGLL[ii,jj] = int(newGridID)
                     gg += 1
                            
       return edgeNodeMapGLL, varCoordGLL, varConGLL

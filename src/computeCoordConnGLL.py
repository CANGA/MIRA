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
    edgePlane *= 1.0 / np.linalg.norm(edgePlane)

    return edgePlane, edgeAngle, RE


def findCoordAlongEdge(x, edgePlane, x1Coord, locAngle, RE):

    # Equation 1 = the angle between x and x1Coord is locAngle
    # Equation 2 = the new grid must be on the same plane as the edge
    # Equation 3 = the point is on a sphere of radius RE
    
    X = np.linalg.norm(x)
    X1 = np.linalg.norm(x1Coord)
    
    tvec = np.cross(x1Coord, x)

    return [np.dot(x, x1Coord) - X * X1 * mt.cos(locAngle),
            np.linalg.norm(np.cross(tvec, edgePlane)),
            x[0]**2 + x[1]**2 + x[2]**2 - RE]


def computeCoordConnGLL(NEL, NGED, NGEL, varCoord, varCon,
                        edgeNodeMap, edgeNodeKDTree, seOrder):

    # NEL = total number of elements
    # NGED = number of grids per edge (new mesh)
    # NGEL = number of grids per element (new mesh)
    # NNG = total number of new grids (besides corners)

    # Get the number of grids in original mesh
    NG = varCoord.shape[1]
    # Get the number of grids/edges per element in original mesh
    NEEL = varCon.shape[1]
    # Get the new number of nodes per element on edges only
    #NGEO = NEEL * (NGED - 1)

    # """ 2nd order method
    if seOrder == 2:
        NGEO = 8
        GN = [-1.0,
              0.0,
              +1.0]

    # """ 4th order method
    if seOrder == 4:
        NGEO = 12
        GN = [-1.0,
              -mt.sqrt(0.2),
              +mt.sqrt(0.2),
              +1.0]

    # Change GN to [0.0 1.0]
    GN = 0.5 * np.add(GN, 1.0)

    # Initialize new global GLL connectivity
    varConGLL = np.zeros((NEL, NGEL))
    # Initialize new global GLL complement of grids
    varCoordGLL = varCoord
    # Initialize new global GLL edge-node map (last column is the cell id)
    NED = edgeNodeMap.shape[0]
    edgeNodeMapGLL = np.zeros((NED, NGED))

    # Loop over the edges
    newGridID = NG + 1
    for ii in range(NED):
        # Get the local node pair map for these edges
        thisEdge = edgeNodeMap[ii, 0:2]

        # Get the opposing edge from neighbor element
        thatEdge = thisEdge[::-1]

        # Set the end points to the new GLL connectivity
        edgeNodeMapGLL[ii, 0] = int(thisEdge[0])
        edgeNodeMapGLL[ii, NGED - 1] = int(thisEdge[1])

        # Check for degenerate edge leaves a 0 in the connectivity
        if thisEdge[0] == thisEdge[1]:
            continue

        # Check that edgeNodeMapGLL hasn't been already set for this edge
        if sum(edgeNodeMapGLL[ii, 1:NGED - 1]) != 0.0:
            continue

        # Fetch the end point coordinates
        coord1 = varCoord[:, int(thisEdge[0]) - 1]
        coord2 = varCoord[:, int(thisEdge[1]) - 1]

        edgePlane, edgeAngle, RE = computeEdgeParametersGLL(
            coord1, coord2)

        # Loop over the new grids on the edge
        for jj in range(1, NGED - 1):
            
            # Compute angular location of new grid on the edge
            newGridAngle = GN[jj] * edgeAngle
            
            # Set the initial guess by rotation about edgePlane unit vector
            vec12 = np.cross(edgePlane, coord1)
            guessCoord = np.dot(edgePlane, coord1) * edgePlane + \
                         mt.cos(newGridAngle) * np.cross(vec12, edgePlane) + \
                         mt.sin(newGridAngle) * vec12
            
            # Solve for the new coordinate
            sol = optimize.root(findCoordAlongEdge, guessCoord,
                (edgePlane, coord1, newGridAngle, RE))
            
            # Store the new grid
            newGrid = np.zeros((3, 1))
            newGrid[:, 0] = sol.x
            varCoordGLL = np.append(varCoordGLL, newGrid, axis=1)
            # Update the new connectivity
            edgeNodeMapGLL[ii, jj] = int(newGridID)
            # Make a new grid ID
            newGridID += 1

        # Set connectivity for the opposing edge in edgeNodeMap (thatEdge)
        edex = edgeNodeKDTree.query_ball_point(
            thatEdge, COINCIDENT_TOLERANCE, p=2, eps=0)
        edgeNodeMapGLL[edex, 1:NGED - 1] = edgeNodeMapGLL[ii, NGED - 2:0:-1]

    # Take edges 2 and 4 and use the interior nodes to build interior of element
    edex = range(NEEL)
    for ii in range(NEL):
        # Fectch every NEEL edges per element
        if ii > 0:
            edex = np.add(edex, NEEL)
        thisElement = edgeNodeMapGLL[edex,:]

        # Set perimeter edges to the new element connectivity
        cdex = range(1, NGED)
        for jj in range(NEEL):
            if jj == 0:
                # First edge in the new connectivity
                varConGLL[ii, 0:NGED] = thisElement[0,:]
            elif jj == NEEL - 1:
                # Last edge in the new connectivity
                if seOrder == 2:
                       varConGLL[ii, (NGEO - 2):NGEO] = thisElement[jj, 0:NGED - 1]
                elif seOrder == 4:
                       varConGLL[ii, (NGEO - 2):NGEO] = thisElement[jj, 1:NGED - 1]
            else:
                cdex = np.add(cdex, (NGED - 1))
                varConGLL[ii, cdex] = thisElement[jj, 1:NGED]

        # Order follows connectivity of the parent element
        if seOrder == 2:
               # Get the ray between the interior node of edges 2 and 4
               edge24 = np.array([[thisElement[3, 1], thisElement[1, 1]]]) - 1
               
               # Fetch the end point coordinates as (2, 3, 2) array
               coords = []
               coords.append((varCoordGLL[:, int(edge24[0, 0])],
                           varCoordGLL[:, int(edge24[0, 1])]))
                         
        elif seOrder == 4:
               # Get the rays between the interior nodes of edges 2 and 4
               edge24 = np.array([[thisElement[3, 2], thisElement[1, 1]],
                                  [thisElement[1, 2], thisElement[3, 1]]]) - 1
               
               # Fetch the end point coordinates as (2, 3, 2) array
               coords = []
               coords.append((varCoordGLL[:, int(edge24[0, 0])],
                           varCoordGLL[:, int(edge24[0, 1])]))
                         
               coords.append((varCoordGLL[:, int(edge24[1, 0])],
                           varCoordGLL[:, int(edge24[1, 1])]))
        else:
               print('ERROR: SE orders other than 2 or 4 NOT supported.')
               return
            
        # Loop over the new interior edges
        for jj in range(len(coords)):
            coord1 = coords[jj][0]
            coord2 = coords[jj][1]
            edgePlane, edgeAngle, RE = computeEdgeParametersGLL(
                coord1, coord2)

            # Loop over the new grids on the edge
            if seOrder == 2:
                   hh = jj
            elif seOrder == 4:
                   hh = 2 * jj
                   
            for kk in range(1, NGED - 1):
                # Compute angular location of new grid on the edge
                newGridAngle = GN[kk] * edgeAngle
                
                # Set the initial guess by rotation about edgePlane unit vector
                vec12 = np.cross(edgePlane, coord1)
                guessCoord = np.dot(edgePlane, coord1) * edgePlane + \
                             mt.cos(newGridAngle) * np.cross(vec12, edgePlane) + \
                             mt.sin(newGridAngle) * vec12
                ''' DEBUG STATEMENT             
                if ii == 0:
                    print(coords[jj][0], newGridAngle, GN[jj], edgeAngle)
                    print(guessCoord)
                '''
                # Solve for the new coordinate
                sol = optimize.root(
                    findCoordAlongEdge, guessCoord, (edgePlane, coord1, newGridAngle, RE))
                
                # Store the new grid
                newGrid = np.zeros((3, 1))
                newGrid[:, 0] = sol.x
                varCoordGLL = np.append(varCoordGLL, newGrid, axis=1)
                # Put the new grid ID at the end of the connectivity row
                varConGLL[ii, NGEO + hh] = newGridID
                # Make a new grid ID (continuing from above)
                newGridID += 1
                hh += 1
    '''
    # VISUAL CHECK OF A GIVEN ELEMENT
    edex = 0
    thisElConn = varConGLL[edex,:].astype(np.int32)
    print('ELEMENT COORDINATES')
    print(varCoordGLL[:,thisElConn-1].T)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for node in thisElConn:
        ndex = node - 1
        point = varCoordGLL[:,ndex]
    
        ax.scatter(point[0], point[1], point[2], marker='D', label=str(node))
        ax.text(point[0], point[1], point[2], str(node))
    plt.legend()
    input('GLL spot check of one element...')
    '''
    return edgeNodeMapGLL, varCoordGLL, varConGLL

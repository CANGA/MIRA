#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:29:22 2018

Computes the area of a given finite volume element of type SHELL4. Uses
connectivity array to index 4 coordinates and reconstruct the area patch

coords: double node data from .g file
connect: int connectivity data

@author: jeguerra
"""
#import math as mt
import numpy as np


def computeAreaWeightPlanarTriangles(coords, connect):

    # Initialize the area
    dFaceArea = 0.0

    # Change connect to 0-based
    cdex = connect - 1
    # Get the coordinates of the quad patch
    nodes = coords[:, cdex]

    # Set the number of subtriangles
    NST = np.size(nodes, axis=1) - 2

    # Loop over the subtriangles and add up the areas
    for ii in range(NST):
        # Gather the coordinate components
        node1 = nodes[:, 0]
        node2 = nodes[:, ii + 1]
        node3 = nodes[:, ii + 2]

        # Compute vectors spanning the triangle patch
        nvec1 = np.subtract(node1, node2)
        nvec2 = np.subtract(node1, node3)

        # Compute the cross product
        ncross = np.cross(nvec1, nvec2)

        dFaceArea += 0.5 * abs(np.linalg.norm(ncross))

    return dFaceArea

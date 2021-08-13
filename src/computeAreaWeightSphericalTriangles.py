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
import math as mt
import numpy as np

def computeAreaWeightSphericalTriangles(coords, connect):
       
       # center node is the first argument
       def computeAngle(cnode, anode, bnode):
              
              # Compute the intersecting radiant plane normals
              cross1 = np.cross(cnode, anode)
              cross2 = np.cross(cnode, bnode)
              
              mc1 = np.linalg.norm(cross1)
              mc2 = np.linalg.norm(cross2)
              
              if (mc1 <= 1.0E-14): 
                     uv1 = 0.0 * cross1
              else:
                     uv1 = 1.0 / mc1 * cross1
                     
              if (mc2 <= 1.0E-14):
                     uv2 = 0.0 * cross2
              else:
                     uv2 = 1.0 / mc2 * cross2
              
              # Compute angle between radiant planes
              uvdot = np.dot(uv1, uv2)
              
              # Fix any round off problems
              if uvdot > 1.0:
                     uvdot = 1.0
              elif uvdot < -1.0:
                     uvdot = -1.0
                     
              angle = mt.acos(uvdot)
              
              return angle
       
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
              node2 = nodes[:, ii+1]
              node3 = nodes[:, ii+2]
              
              # Compute 3 included angles
              angle1 = computeAngle(node1, node2, node3)
              angle2 = computeAngle(node2, node1, node3)
              angle3 = computeAngle(node3, node1, node2)
              
              # Compute the area by Girard's Theorem
              dFaceArea += (angle1 + angle2 + angle3 - mt.pi)
                            
       return dFaceArea
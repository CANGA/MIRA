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
from scipy.linalg import norm

# Order 4 Gauss quadrature nodes and weights
def getGaussNodesWeights(order):
       #""" 2nd order method
       if order == 2:
              GN = [-1.0, \
                     0.0
                    +1.0]
              
              GW = [1.0 / 3.0, \
                    4.0 / 3.0, \
                    1.0 / 3.0]
       #""" 4th order method
       if order == 4:
              GN = [-1.0, \
                    -mt.sqrt(0.2), \
                    +mt.sqrt(0.2), \
                    +1.0]
              
              GW = [1.0 / 6.0, \
                    5.0 / 6.0, \
                    5.0 / 6.0, \
                    1.0 / 6.0]
       #"""
       
       return np.ravel(GN), \
              np.ravel(GW)
              
def computeCart2LL(cellCoord):

       RO = np.linalg.norm(cellCoord)
       psi = mt.asin(1.0 / RO * cellCoord[2])
       lam = mt.atan2(-cellCoord[0], -cellCoord[1]) + mt.pi
       pointLonLat = [360.0 * lam / (2.0 * mt.pi), 180.0 * psi / mt.pi]
              
       return pointLonLat

def computeSEIntegral(varGLL, nodes, order, avg, farea):
       # avg = Boolean flag to take average of the function
       # farea = Boolean flag to compute only the area integral (ignore field)
       
       # Initialize the area
       dFaceArea = 0.0
       
       # Initialize the integral
       dFunIntegral = 0.0
       
       # Loop over the subtriangles and add up the areas
       GN, GW = getGaussNodesWeights(order)
       NP = len(GW)
       
       # nodes is a 3D array [3, NP, NP] of the element's GLL grid (in global coords)

       # Loop over the quadrature points
       for pp in range(NP):
              for qq in range(NP):

                     # Fetch global coords of this quadrature point
                     dF = nodes[:,pp,qq]
                     dF2 = dF**2
                     dR = norm(dF, 2)

                     # Sample SH field at this quadrature point
                     # Convert dF to Lon/Lat
                     if farea == False:
                            thisVar = varGLL
                     elif farea == True:
                            thisVar = 1.0
                     
                     # Local difference between this node and the next in a direction
                     dDaF = dnDiffa
                     # Local difference between this node and the next in b direction
                     dDbF = dnDiffb

                     dR = norm(dF, 2)

                     dDGMat = np.array([[(dF2[1] + dF2[2]), - dF[0] * dF[1], - dF[0] * dF[2]], \
                                        [- dF[1] * dF[0], (dF2[0] + dF2[2]), - dF[1] * dF[2]], \
                                        [- dF[2] * dF[0], - dF[2] * dF[1], (dF2[0] + dF2[1])]]) 
                     
                     dDaG = np.dot(dDGMat, dDaF)
                     dDbG = np.dot(dDGMat, dDbF)

                     dDenomTerm = 1.0 / (dR**3)
                     
                     # This happens to dDaG twice...
                     dDaG *= dDenomTerm
                     dDbG *= dDenomTerm
                     
                     dJV = np.cross(dDaG, dDbG)
                     dJacobianGWppqq = norm(dJV, 2) * GW[pp] * GW[qq]
                     
                     # Sum up the cell area
                     dFaceArea += dJacobianGWppqq
                     # Sum up the integral of the field
                     dFunIntegral += thisVar * dJacobianGWppqq
       
       
       # Compute the cell average                            
       dFunAverage = dFunIntegral / dFaceArea
       
       # When a cell average is required
       if avg:
              dFunIntegral = dFunAverage
                            
       return dFunIntegral


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
       
       # Scale the points/weights to [0 1]
       ovec = np.ones(np.size(GN))
       GN = 0.5 * np.matrix(np.add(GN, ovec))
       GW = 0.5 * np.matrix(GW)
       
       return np.ravel(GN), \
              np.ravel(GW)
              
def computeCart2LL(cellCoord):

       RO = np.linalg.norm(cellCoord)
       psi = mt.asin(1.0 / RO * cellCoord[2])
       lam = mt.atan2(-cellCoord[0], -cellCoord[1]) + mt.pi
       pointLonLat = [360.0 * lam / (2.0 * mt.pi), 180.0 * psi / mt.pi]
              
       return pointLonLat

def computeSEIntegral(varGLL, varCoords, order, avg, farea):
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
       nodes = np.zeros((3,NP,NP))
       cdex = [0, 1, 2, 3, 11, 12, 13, 4, 10, 15, 14, 5, 9, 8, 7, 6]
       rr = 0
       for pp in range(NP):
              for qq in range(NP):
                     nodes[:,pp,qq] = varCoords[:,cdex[rr]]
                     rr += 1
       
       # Compute the plane edge directions for this quadrilateral
       nD30 = np.subtract(varCoords[:,3], varCoords[:,0])
       nD69 = np.subtract(varCoords[:,6], varCoords[:,9])
       nD90 = np.subtract(varCoords[:,9], varCoords[:,0])
       nD63 = np.subtract(varCoords[:,6], varCoords[:,3])

       # Loop over the quadrature points
       dF = np.zeros((3,))
       dDaF = np.zeros((3,1))
       dDbF = np.zeros((3,1))
       for pp in range(NP):
              for qq in range(NP):
                     # Enhance...
                     dA = GN[pp]
                     dB = GN[qq]
                     
                     dOmdA = (1.0 - dA)
                     dOmdB = (1.0 - dB)

                     # Fetch global coords of this quadrature point (on the plane)
                     dF[0] = varCoords[0,0] * dOmdA * dOmdB \
              		+ varCoords[0,3] * dA * dOmdB \
              		+ varCoords[0,6] * dA * dB \
              		+ varCoords[0,9] * dOmdA * dB
                     dF[1] = varCoords[1,0] * dOmdA * dOmdB \
              		+ varCoords[1,3] * dA * dOmdB \
              		+ varCoords[1,6] * dA * dB \
              		+ varCoords[1,9] * dOmdA * dB
                     dF[2] = varCoords[2,0] * dOmdA * dOmdB \
              		+ varCoords[2,3] * dA * dOmdB \
              		+ varCoords[2,6] * dA * dB \
              		+ varCoords[2,9] * dOmdA * dB
                            
                     dF2 = dF**2
                     dR = norm(dF, 2)

                     # Sample SH field at this quadrature point
                     # Convert dF to Lon/Lat
                     if farea == False:
                            thisVar = varGLL[pp,qq]
                     elif farea == True:
                            thisVar = 1.0
                     
                     # Local secant vector between this node and the next in a direction
                     dDaF[:,0] = [dOmdB * nD30[0] + dB * nD69[0], \
                             dOmdB * nD30[1] + dB * nD69[1], \
                             dOmdB * nD30[2] + dB * nD69[2]]
                     
                     # Local secant vector between this node and the next in b direction
                     dDbF[:,0] = [dOmdA * nD90[0] + dA * nD63[0], \
                             dOmdA * nD90[1] + dA * nD63[1], \
                             dOmdA * nD90[2] + dA * nD63[2]]
                     
                     # Set the spherical metric...
                     dDGMat = np.array([[(dF2[1] + dF2[2]), - dF[0] * dF[1], - dF[0] * dF[2]], \
                                        [- dF[1] * dF[0], (dF2[0] + dF2[2]), - dF[1] * dF[2]], \
                                        [- dF[2] * dF[0], - dF[2] * dF[1], (dF2[0] + dF2[1])]]) 
                     
                     # Compute the local tangent vectors in a and b directions
                     dDaG = np.dot(dDGMat, dDaF)
                     dDbG = np.dot(dDGMat, dDbF)

                     dDenomTerm = 1.0 / (dR**3)
                     
                     # This happens to dDaG twice...
                     dDaG *= dDenomTerm
                     dDbG *= dDenomTerm
                     
                     dJV = np.cross(np.ravel(dDaG), np.ravel(dDbG))
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


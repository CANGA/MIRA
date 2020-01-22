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
       """
       # 2nd order method
       if order == 2:
              GN = [-0.5773502691896257, \
                    +0.5773502691896257]
              
              GW = [+1.0, \
                    +1.0]
       # 4th order method
       if order == 4:
              GN = [-0.8611363115940526, \
                    -0.3399810435848563, \
                    +0.3399810435848563, \
                    +0.8611363115940526]
              
              GW = [0.3478548451374538, \
                    0.6521451548625461, \
                    0.6521451548625461, \
                    0.3478548451374538]
       #
       # 6th order method
       if order == 6:
              GN = [-0.9324695142031521, \
                    -0.6612093864662645, \
                    -0.2386191860831969, \
                    +0.2386191860831969, \
                    +0.6612093864662645, \
                    +0.9324695142031521]
              
              GW = [0.1713244923791704, \
                    0.3607615730481386, \
                    0.4679139345726910, \
                    0.4679139345726910, \
                    0.3607615730481386, \
                    0.1713244923791704]
       #
       
       ###
       """

       from numpy.polynomial.legendre import leggauss
       GN, GW = leggauss(order)

       # Scale the points/weights to [0 1]
       ovec = np.ones(np.size(GN))
       GN = 0.5 * np.matrix(np.add(GN, ovec))
       GW = 0.5 * np.matrix(GW)
      
       # print('Quadratures: ', GN, GW)
       return np.ravel(GN), \
              np.ravel(GW)
              
def computeCart2LL(cellCoord):

       RO = np.linalg.norm(cellCoord)
       psi = mt.asin(1.0 / RO * cellCoord[2])
       lam = mt.atan2(-cellCoord[0], -cellCoord[1]) + mt.pi
       pointLonLat = [360.0 * lam / (2.0 * mt.pi), 180.0 * psi / mt.pi]
              
       return pointLonLat

def computeAreaIntegral(clm, nodes, order, avg, farea):
       # avg = Boolean flag to take average of the function
       # farea = Boolean flag to compute only the area integral (ignore field)
       
       # Initialize the area
       dFaceArea = 0.0
       
       # Initialize the integral
       dFunIntegral = 0.0
       
       # Set the number of subtriangles
       NST = np.size(nodes, axis=1) - 2
       
       # Loop over the subtriangles and add up the areas
       GN, GW = getGaussNodesWeights(order)
       NP = len(GW)

       for ii in range(NST):
              # Gather the coordinate components
              node1 = nodes[:,0]
              node2 = nodes[:,ii+1]
              node3 = nodes[:,ii+2]

              nDMatrix = np.array([ [node1[0], node1[1], node1[2]],
                                     [node2[0], node2[1], node2[2]],
                                     [node3[0], node3[1], node3[2]] ])

              nD21 = np.array([ node2[0]-node1[0], node2[1]-node1[1], node2[2]-node1[2] ])

              # Loop over the quadrature points
              for pp in range(NP):
                     for qq in range(NP):
                            # Enhance...
                            dA = GN[pp]
                            dB = GN[qq]
                            
                            dOmdA = (1.0 - dA)
                            dOmdB = (1.0 - dB)

                            dAOmdA = np.array([ [-dOmdA], [-dA], [1.0] ])
                            dOmdBOmdA = np.array([ dOmdB * dOmdA, dOmdB * dA, dB ])

                            # Compute global coords of this quadrature point
                            dF = np.dot(nDMatrix.T, dOmdBOmdA)
                            dF2 = dF**2

                            # Sample SH field at this quadrature point
                            # Convert dF to Lon/Lat
                            if farea == False:
                                   dFLonLat = computeCart2LL(dF)
                                   thisVar = clm.expand(lon=dFLonLat[0], lat=dFLonLat[1])
                            elif farea == True:
                                   thisVar = 1.0

                            dDaF = (dOmdB * nD21)

                            dDbF = (np.dot(nDMatrix.T, dAOmdA).T)[0]

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


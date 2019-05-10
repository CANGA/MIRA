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

# Order 4 Gauss quadrature nodes and weights
def getGaussNodesWeights(order):
       if order == 2:
              GN = [-0.5773502691896257, \
                    +0.5773502691896257]
              
              GW = [+1.0, \
                    +1.0]
       #""" 4th oder method for testing
       if order == 4:
              GN = [-0.8611363115940526, \
                    -0.3399810435848563, \
                    +0.3399810435848563, \
                    +0.8611363115940526]
              
              GW = [0.3478548451374538, \
                    0.6521451548625461, \
                    0.6521451548625461, \
                    0.3478548451374538]
       #"""
       #""" 6th order method is slower
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
              
              n1x = node1[0]
              n1y = node1[1]
              n1z = node1[2]
              n2x = node2[0]
              n2y = node2[1]
              n2z = node2[2]
              n3x = node3[0]
              n3y = node3[1]
              n3z = node3[2]
              
              # Loop over the quadrature points
              for pp in range(NP):
                     for qq in range(NP):
                            # Enhance...
                            dA = GN[pp]
                            dB = GN[qq]
                            
                            dOmdA = (1.0 - dA)
                            dOmdB = (1.0 - dB)
                            
                            # Compute global coords of this quadrature point
                            dF = [dOmdB * (dOmdA * n1x + dA * n2x) + dB * n3x, \
                                  dOmdB * (dOmdA * n1y + dA * n2y) + dB * n3y, \
                                  dOmdB * (dOmdA * n1z + dA * n2z) + dB * n3z]
                            
                            # Sample SH field at this quadrature point
                            # Convert dF to Lon/Lat
                            if farea == False:
                                   dFLonLat = computeCart2LL(dF)
                                   thisVar = clm.expand(lon=dFLonLat[0], lat=dFLonLat[1])
                            elif farea == True:
                                   thisVar = 1.0
                            
                            dDaF = [dOmdB * (n2x - n1x), \
                                    dOmdB * (n2y - n1y), \
                                    dOmdB * (n2z - n1z)]
                            
                            dDbF = [-dOmdA * n1x - dA * n2x + n3x, \
                                    -dOmdA * n1y - dA * n2y + n3y, \
                                    -dOmdA * n1z - dA * n2z + n3z]
                            
                            dR = mt.sqrt(dF[0]**2 + dF[1]**2 + dF[2]**2)
                            
                            dDaG = [dDaF[0] * (dF[1]**2 + dF[2]**2) - dF[0] * (dDaF[1] * dF[1] + dDaF[2] * dF[2]), \
                                    dDaF[1] * (dF[0]**2 + dF[2]**2) - dF[1] * (dDaF[0] * dF[0] + dDaF[2] * dF[2]), \
                                    dDaF[2] * (dF[0]**2 + dF[1]**2) - dF[2] * (dDaF[0] * dF[0] + dDaF[1] * dF[1])]
                            
                            dDbG = [dDbF[0] * (dF[1]**2 + dF[2]**2) - dF[0] * (dDbF[1] * dF[1] + dDbF[2] * dF[2]), \
                                    dDbF[1] * (dF[0]**2 + dF[2]**2) - dF[1] * (dDbF[0] * dF[0] + dDbF[2] * dF[2]), \
                                    dDbF[2] * (dF[0]**2 + dF[1]**2) - dF[2] * (dDbF[0] * dF[0] + dDbF[1] * dF[1])]
                            
                            dDenomTerm = 1.0 / (dR**3)
                            
                            # This happens to dDaG twice...
                            dDaG = dDenomTerm * np.matrix(dDaG)
                            dDbG = dDenomTerm * np.matrix(dDbG)
                            
                            dJV = np.cross(dDaG, dDbG)
                            dJacobian = mt.sqrt(dJV[0,0]**2 + \
                                                dJV[0,1]**2 + \
                                                dJV[0,2]**2)
                            
                            # Sum up the cell area
                            dFaceArea += GW[pp] * GW[qq] * dJacobian
                            # Sum up the integral of the field
                            dFunIntegral += thisVar * GW[pp] * GW[qq] * dJacobian
       
       
       # Compute the cell average                            
       dFunAverage = dFunIntegral / dFaceArea
       
       # When a cell average is required
       if avg:
              dFunIntegral = dFunAverage
                            
       return dFunIntegral
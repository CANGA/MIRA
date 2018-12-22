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

def computeAreaWeight(coords, connect):
       
       # Order 4 Gauss quadrature nodes and weights
       # TO DO: Add order 6 as done in tempest remap
       def getGaussNodesWeights():
              GN = [-0.8611363115940526, \
                    -0.3399810435848563, \
                    +0.3399810435848563, \
                    +0.8611363115940526]
              
              GW = [0.3478548451374538, \
                    0.6521451548625461, \
                    0.6521451548625461, \
                    0.3478548451374538]
              
              # Scale the points/weights to [0 1]
              GN = 1.0 + 0.5 * np.matrix(GN)
              GW = 0.5 * np.matrix(GW)
              
              return np.ravel(GN), \
                     np.ravel(GW)
       
       # Initialize the area
       dFaceArea = 0.0
       
       # Change connect to 0-based
       cdex = connect - 1
       # Get the coordinates of the quad patch
       nodes = coords[:, cdex]
       
       # Set the number of subtriangles
       NST = len(nodes) - 2
       
       # Loop over the subtriangles and add up the areas
       NP = 4
       GN, GW = getGaussNodesWeights()
       for ii in range(NST):
              # Gather the coordinate components
              node1 = nodes[:, 0]
              node2 = nodes[:, ii+1]
              node3 = nodes[:, ii+2]
              
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
                            
                            dOmdB = (1.0 - dB)
                            dOmdA = (1.0 - dA)
                            
                            dF = [dOmdB * (dOmdA * n1x + dA * n2x) + dB * n3x, \
                                  dOmdB * (dOmdA * n1y + dA * n2y) + dB * n3y, \
                                  dOmdB * (dOmdA * n1z + dA * n2z) + dB * n3z]
                            
                            dDaF = [dOmdB * (n1x - n2x), \
                                    dOmdB * (n2y - n1y), \
                                    dOmdB * (n2z - n1z)]
                            
                            dDbF = [-dOmdA * n1x - dA * n2x  + n3x, \
                                    -dOmdA * n1y - dA * n2y  + n3y, \
                                    -dOmdA * n1z - dA * n2z  + n3z]
                            
                            dR = mt.sqrt(dF[0]**2 + dF[1]**2 + dF[2]**2)
                            
                            dDaG = [dDaF[0] * (dF[1]**2 + dF[2]**2) - dF[0] * (dDaF[1] * dF[1] + dDaF[2] * dF[2]), \
                                    dDaF[1] * (dF[0]**2 + dF[2]**2) - dF[1] * (dDaF[0] * dF[0] + dDaF[2] * dF[2]), \
                                    dDaF[2] * (dF[0]**2 + dF[1]**2) - dF[2] * (dDaF[0] * dF[0] + dDaF[1] * dF[1])]
                            
                            dDbG = [dDbF[0] * (dF[1]**2 + dF[2]**2) - dF[0] * (dDbF[1] * dF[1] + dDbF[2] * dF[2]), \
                                    dDbF[1] * (dF[0]**2 + dF[2]**2) - dF[1] * (dDbF[0] * dF[0] + dDbF[2] * dF[2]), \
                                    dDbF[2] * (dF[0]**2 + dF[1]**2) - dF[2] * (dDbF[0] * dF[0] + dDbF[1] * dF[1])]
                            
                            dDenomTerm = 1.0 / (dR**3)
                            
                            dDaG = dDenomTerm * np.matrix(dDaG)
                            dDbG = dDenomTerm * np.matrix(dDbG)
                            
                            dJV = np.cross(dDaG, dDbG)
                            dJacobian = mt.sqrt(dJV[0,0]**2 + \
                                                dJV[0,1]**2 + \
                                                dJV[0,2]**2)
                            
                            dFaceArea += GW[pp] * GW[qq] * dJacobian
                            
       return dFaceArea
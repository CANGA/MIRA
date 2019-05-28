#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 13:35:59 2019

@author: jeguerra
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 24 08:26:47 2018

Compute a gradient on the GLL grid.

By construction, all the vectors involved start at the origin... convenient.
Parameterized flux integral assumes constant radius and field value.

@author: jeguerra
"""

import numpy as np
import math as mt

# Order 4 Gauss quadrature nodes
def getGLLNodesWeights(order):
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

def computeTangentBasisSE(varCoords, order):
       
       # Loop over the subtriangles and add up the areas
       GN, GW = getGLLNodesWeights(order)
       NP = len(GW)
       
       # Set the connectivity index vector corresponding to varCoords
       cdex = [0, 1, 2, 3, 11, 12, 13, 4, 10, 15, 14, 5, 9, 8, 7, 6]
       
       # Compute the plane edge directions for this quadrilateral
       nD30 = np.subtract(varCoords[:,3], varCoords[:,0])
       nD69 = np.subtract(varCoords[:,6], varCoords[:,9])
       nD90 = np.subtract(varCoords[:,9], varCoords[:,0])
       nD63 = np.subtract(varCoords[:,6], varCoords[:,3])

       # Loop over the quadrature points
       dF = np.zeros((3,))
       dDaF = np.zeros((3,1))
       dDbF = np.zeros((3,1))
       dDaG = np.zeros(3,(len(cdex)))
       dDbG = np.zeros(3,(len(cdex)))
       rr = 0
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
                     dR = np.linalg.norm(dF, 2)
                     
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
                     dDaG[cdex[rr]] = np.dot(dDGMat, dDaF)
                     dDbG[cdex[rr]] = np.dot(dDGMat, dDbF)

                     dDenomTerm = 1.0 / (dR**3)
                     
                     # This happens to dDaG twice...
                     dDaG[cdex[rr]] *= dDenomTerm
                     dDbG[cdex[rr]] *= dDenomTerm
                     rr += 1
                     
       # Return the local tangent basis vectors
       return dDaG, dDbG

def computeDerivativeMatrixGLL(order):
       
       if order == 4:
              GN, GW = getGLLNodesWeights(order)
              
              # Derivatives of Legendre polynomials to order 4
              DD = np.array([[0.0, 1.0, 3.0 * GN[0], 0.5 * (15.0 * GN[0]** - 3.0)], \
                             [0.0, 1.0, 3.0 * GN[1], 0.5 * (15.0 * GN[1]** - 3.0)], \
                             [0.0, 1.0, 3.0 * GN[2], 0.5 * (15.0 * GN[2]** - 3.0)], \
                             [0.0, 1.0, 3.0 * GN[3], 0.5 * (15.0 * GN[3]** - 3.0)]])
              
              # Scale derivative matrix to [0.0 1.0]
              DD *= 2.0
              
       return DD

def computeGradientSE(varList, varCon, varCoords, order, jacobians):
       SF = np.float64
       
       # Gradients are 3 component vectors
       nc = 3
       gradShape1 = (nc, varList[0].shape[0])
       gradShape2 = (nc, varList[1].shape[0])
       
       # Set the main return variable
       varGradient = [np.zeros(gradShape1, dtype=SF), \
                      np.zeros(gradShape2, dtype=SF)]
       
       # Alpha direction connectivity index (columns)
       cdex = [0, 1, 2, 3, 11, 12, 13, 4, 10, 15, 14, 5, 9, 8, 7, 6]
       # Beta direction connectivity index (rows)
       rdex = [0, 11, 10, 9, 1, 12, 15, 8, 2, 13, 14, 7, 3, 4, 5, 6]
       
       NC = int(varCon.shape[0]) # number of elements
       NP = order
       
       # Get the reference domain derivative matrix at the given order
       dDab = computeDerivativeMatrixGLL(order)
       
       # Compute the local tangent basis vectors
       dDaG, dDbG = computeTangentBasisSE(varCoords, order)
       
       varST = np.zeros((order,1))
       varS2T = np.zeros((order,1))
       varDervAlphaST = np.zeros(len(cdex))
       varDervAlphaS2T = np.zeros(len(cdex))
       varDervBetaST = np.zeros(len(rdex))
       varDervBetaS2T = np.zeros(len(rdex))
       # Loop over the elements
       for jj in range(NC):
              
              ndex = varCon[jj,:].astype(int) - 1
              
              # Fetch the variable on this element
              var1 = varList[0][ndex]
              var2 = varList[1][ndex]
              
              # Loop over GLL order and compute local derivatives
              for oo in range(NP):
                     # Get the alpha index
                     odex = cdex[0:order] + NP * oo
                     # Get variables as an (order X 1) vector
                     varST[:,0] = var1[odex.astype(int)]
                     varS2T[:,0] = var2[odex.astype(int)]
                     # Compute the native derivative in alpha
                     varDervAlphaST[odex.astype(int)] = np.dot(dDab, varST)
                     varDervAlphaS2T[odex.astype(int)] = np.dot(dDab, varS2T)
                     
                     # Get the beta index
                     odex = rdex[0:order] + NP * oo
                     # Get variables as an (order X 1) vector
                     varST[:,0] = var1[odex.astype(int)]
                     varS2T[:,0] = var2[odex.astype(int)]
                     # Compute the native derivative in beta
                     varDervBetaST[odex.astype(int)] = np.dot(dDab, varST)
                     varDervBetaS2T[odex.astype(int)] = np.dot(dDab, varS2T)
                     
              # Set the gradient to the appropriate grids
              varGradient[0][:,ndex] = np.mul(varDervAlphaST, dDaG, axis=1)
              varGradient[0][:,ndex] += np.mul(varDervBetaST, dDbG, axis=1)
              
              varGradient[1][:,ndex] = np.mul(varDervAlphaS2T, dDaG, axis=1)
              varGradient[1][:,ndex] += np.mul(varDervBetaS2T, dDbG, axis=1)
              
              
       return varGradient
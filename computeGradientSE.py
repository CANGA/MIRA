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
import sys
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
       dDaG = np.zeros((3,len(cdex)))
       dDbG = np.zeros((3,len(cdex)))
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
                     dDaG[:,cdex[rr]] = np.ravel(np.dot(dDGMat, dDaF))
                     dDbG[:,cdex[rr]] = np.ravel(np.dot(dDGMat, dDbF))

                     dDenomTerm = 1.0 / (dR**3)
                     
                     # This happens to dDaG twice...
                     dDaG[:,cdex[rr]] *= dDenomTerm
                     dDbG[:,cdex[rr]] *= dDenomTerm
                     rr += 1
                     
       # Return the local tangent basis vectors
       return dDaG, dDbG

def computeDerivativeMatrixGLL(order):
       
       if order == 4:
              GN, GW = getGLLNodesWeights(order)
              
              # Derivatives of Legendre polynomials to order 4
              DD = np.array([[0.0, 1.0, 3.0 * GN[0], 0.5 * (15.0 * GN[0]**2.0 - 3.0)], \
                             [0.0, 1.0, 3.0 * GN[1], 0.5 * (15.0 * GN[1]**2.0 - 3.0)], \
                             [0.0, 1.0, 3.0 * GN[2], 0.5 * (15.0 * GN[2]**2.0 - 3.0)], \
                             [0.0, 1.0, 3.0 * GN[3], 0.5 * (15.0 * GN[3]**2.0 - 3.0)]])
              
              # Scale derivative matrix to [0.0 1.0]
              DD *= 2.0
              
       return DD

def computeGradientSE(varField, varCon, varCoords, order, jacobians):
       
       # Check GLL connectivity
       if int(varCon.shape[1]) != int(jacobians.shape[1]):
              print('ERROR: INVALID MESH DATA: Jacobians DOF mismatch connectivity DOF!')
              sys.exit()
       
       SF = np.float64
       
       # Gradients are 3 component vectors
       nc = 3
       gradShape = (nc, varField.shape[0])
       
       # Set the main return variable
       varGradient = np.zeros(gradShape, dtype=SF)
       
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
       
       varT = np.zeros((order,1))
       varDervAlphaT = np.zeros(len(cdex))
       varDervBetaT = np.zeros(len(rdex))
       
       # Loop over the elements
       for jj in range(NC):
              
              ndex = np.subtract(varCon[jj,:].astype(int), 1)
              
              # Fetch the variable on this element
              varT = varField[ndex]
              
              # Loop over GLL order and compute local derivatives
              for oo in range(NP):
                     
                     # Get the alpha index
                     adex = cdex[(NP * oo):(NP * (oo+1))]
                     # Get variables as an (order X 1) vector
                     varT = varField[adex]
                     # Compute the native derivative in alpha
                     varDervAlphaT[adex] = np.ravel(np.dot(dDab, varT))
                     
                     # Get the beta index
                     bdex = rdex[(NP * oo):(NP * (oo+1))]
                     # Get variables as an (order X 1) vector
                     varT = varField[bdex]
                     # Compute the native derivative in beta
                     varDervBetaT[bdex] = np.ravel(np.dot(dDab, varT))
                     
              # Set the gradient to the appropriate grids             
              varGradient[:,ndex] = np.multiply(varDervAlphaT, dDaG)
              varGradient[:,ndex] += np.multiply(varDervBetaT, dDbG)
                                          
       return varGradient
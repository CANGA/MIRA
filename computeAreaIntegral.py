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
      
       if order < 1:
              print('INVALID QUADRATURE ORDER! SETTING ORDER 2.')
              order = 2
       elif order > 100:
              print('QUADRATURE ORDER GREATER THAN 100... ARE YOU SURE?')

       # Arbitrary order Numpy implementation (reportedly good up to order = 200)
       GN, GW = np.polynomial.legendre.leggauss(order)

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

def compute_dFs(GN, GW, nDMatrix):
    NP = len(GW)
    dOmdBOmdAs = np.zeros(shape=(NP*NP,3), dtype='f8')
    # Loop over the quadrature points
    for pp in range(NP):
           for qq in range(NP):
                  # Enhance...
                  dA = GN[pp]
                  dB = GN[qq]
                  
                  dOmdA = (1.0 - dA)
                  dOmdB = (1.0 - dB)

                  dAOmdA = np.array([ [-dOmdA], [-dA], [1.0] ])
                  dOmdBOmdAs[pp*NP+qq,:] = np.array([ dOmdB * dOmdA, dOmdB * dA, dB ])

    dFs = np.dot(nDMatrix.T, dOmdBOmdAs.T).T
    return dFs

def computeLonLats(GN, GW, dFs):
    NP = len(GW)
    for i in range(NP*NP):
           lonlats[i,:] = computeCart2LL(dFs[i, :])
    return lonlats

def computeAreaIntegral(clm, nodes, GN, GW, avg, farea):
       # avg = Boolean flag to take average of the function
       # farea = Boolean flag to compute only the area integral (ignore field)
       
       # Initialize the area
       dFaceArea = 0.0
       
       # Initialize the integral
       dFunIntegral = 0.0
       
       # Set the number of subtriangles
       NST = np.size(nodes, axis=1) - 2
       
       # Loop over the subtriangles and add up the areas
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

#              dOmdBOmdAs = np.zeros(shape=(NP*NP,3), dtype='f8')
#              dFs = compute_dFs(GN, GW, nDMatrix, dOmdBOmdAs)
#              #lonlats = computeLonLats(GN, GW, dFs)
#              lonlats = np.zeros(shape=(NP*NP,2), dtype='f8')
#              computeLonLats(GN, GW, dFs, lonlats)
#
#              allVar = clm.expand(lon=lonlats[:,0], lat=lonlats[:,1])
#
#              # Loop over the quadrature points
#              for pp in range(NP):
#                     for qq in range(NP):
#                            # Enhance...
#                            dA = GN[pp]
#                            dB = GN[qq]
#                            
#                            dOmdA = (1.0 - dA)
#                            dOmdB = (1.0 - dB)
#
#                            dAOmdA = np.array([ [-dOmdA], [-dA], [1.0] ])
#
#                            # global coords of this quadrature point
#                            dF = dFs[pp*NP+qq,:]
#                            dF2 = dF**2
#
#                            dDaF = (dOmdB * nD21)
#
#                            dDbF = (np.dot(nDMatrix.T, dAOmdA).T)[0]
#
#                            dR = np.linalg.norm(dF, 2)
#
#                            dDGMat = np.array([[(dF2[1] + dF2[2]), - dF[0] * dF[1], - dF[0] * dF[2]], \
#                                               [- dF[1] * dF[0], (dF2[0] + dF2[2]), - dF[1] * dF[2]], \
#                                               [- dF[2] * dF[0], - dF[2] * dF[1], (dF2[0] + dF2[1])]]) 
#                            
#                            dDaG = np.dot(dDGMat, dDaF)
#                            dDbG = np.dot(dDGMat, dDbF)
#
#                            dDenomTerm = 1.0 / (dR**3)
#                            
#                            # This happens to dDaG twice...
#                            dDaG *= dDenomTerm
#                            dDbG *= dDenomTerm
#                            
#                            dJV = np.cross(dDaG, dDbG)
#                            dJacobianGWppqq = np.linalg.norm(dJV, 2) * GW[pp] * GW[qq]
#                            
#                            # Sum up the cell area
#                            dFaceArea += dJacobianGWppqq
#
#                            # Sample SH field at this quadrature point
#                            # Convert dF to Lon/Lat
#                            if farea == False:
#                                   #dFLonLat = computeCart2LL(dF)
#                                   #thisVar = clm.expand(lon=dFLonLat[0], lat=dFLonLat[1])
#                                   thisVar = allVar[pp*NP+qq]
#                            elif farea == True:
#                                   thisVar = 1.0
#
#                            # Sum up the integral of the field
#                            dFunIntegral += thisVar * dJacobianGWppqq
#
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

                            dR = np.linalg.norm(dF, 2)

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
                            dJacobianGWppqq = np.linalg.norm(dJV, 2) * GW[pp] * GW[qq]
                            
                            # Sum up the cell area
                            dFaceArea += dJacobianGWppqq
                            # Sum up the integral of the field
                            dFunIntegral += thisVar * dJacobianGWppqq

       
       
       # When a cell average is required
       if avg:
              # Compute the cell average                            
              dFunAverage = dFunIntegral / dFaceArea
              dFunIntegral = dFunAverage
                            
       return dFunIntegral


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
import computeSphericalCartesianTransforms as trans
import compute_sphere_int as csi
from numba import jit

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


def computeAreaIntegral(clm, nodes, order, avg, farea):

    useGQscheme = True

    if useGQscheme:
        GN, GW = getGaussNodesWeights(order)

        [gq_cell_int, gq_areas] = computeAreaIntegralWithGQ(
            clm, nodes, GN, GW, avg, farea)
        return gq_areas if farea else gq_cell_int

    else:
        surfs = np.array([[0, 1, 2, 3]], dtype=int)
        cell_int, areas = csi.compute_sphere_int(nodes.T, surfs, clm)
        return np.sum(areas) if farea else np.sum(cell_int / areas)


@jit(nopython=True, parallel=False)
def computeAreaIntegralInternal(
        NST, NP, GW, GN, nodes, farea, all_dF, all_dJacobianGWppqq):

    # Initialize the area
    dFaceArea = 0.0

    # Initialize the integral
    dFunIntegral = 0.0

    # Link: https://people.sc.fsu.edu/~jburkardt/f_src/stripack/stripack.f90
    # Reference for function "areas" as an alternate implementation for cell areas
    # on a spherical mesh

    for ii in range(NST):
        # Gather the coordinate components
        node1 = nodes[:, 0]
        node2 = nodes[:, ii + 1]
        node3 = nodes[:, ii + 2]

        nDMatrix = np.array([[node1[0], node2[0], node3[0]],
                             [node1[1], node2[1], node3[1]],
                             [node1[2], node2[2], node3[2]]])

        nD21 = np.array(
            [(node2[0] - node1[0]), (node2[1] - node1[1]), (node2[2] - node1[2])])

        # Loop over the quadrature points
        for pp in range(NP):
            for qq in range(NP):
                # Enhance...
                dA = GN[pp]
                dB = GN[qq]

                dOmdA = (1.0 - dA)
                dOmdB = (1.0 - dB)

                dAOmdA = np.array([-dOmdA, -dA, 1.0])
                dOmdBOmdA = np.array([dOmdB * dOmdA, dOmdB * dA, dB])

                # Compute global coords of this quadrature point
                dF = nDMatrix @ dOmdBOmdA

                dDaF = (dOmdB * nD21)
                dDbF = nDMatrix @ dAOmdA

                dR = np.linalg.norm(dF, 2)

                dDenomTerm = 1.0 / (dR**3)

                dDaG = np.array([
                    dDaF[0] * (dF[1] * dF[1] + dF[2] * dF[2]) -
                    dF[0] * (dDaF[1] * dF[1] + dDaF[2] * dF[2]),
                    dDaF[1] * (dF[0] * dF[0] + dF[2] * dF[2]) -
                    dF[1] * (dDaF[0] * dF[0] + dDaF[2] * dF[2]),
                    dDaF[2] * (dF[0] * dF[0] + dF[1] * dF[1]) -
                    dF[2] * (dDaF[0] * dF[0] + dDaF[1] * dF[1])
                ]) * dDenomTerm

                dDbG = np.array([
                    dDbF[0] * (dF[1] * dF[1] + dF[2] * dF[2]) -
                    dF[0] * (dDbF[1] * dF[1] + dDbF[2] * dF[2]),
                    dDbF[1] * (dF[0] * dF[0] + dF[2] * dF[2]) -
                    dF[1] * (dDbF[0] * dF[0] + dDbF[2] * dF[2]),
                    dDbF[2] * (dF[0] * dF[0] + dF[1] * dF[1]) -
                    dF[2] * (dDbF[0] * dF[0] + dDbF[1] * dF[1])
                ]) * dDenomTerm

                dJV = np.cross(dDaG, dDbG)
                dJacobianGWppqq = np.linalg.norm(dJV, 2) * GW[pp] * GW[qq]

                # Sum up the cell area
                dFaceArea += dJacobianGWppqq

                # Sample SH field at this quadrature point
                if farea:
                    # Sum up the integral of the field
                    dFunIntegral += dJacobianGWppqq
                else:
                    # Convert dF to Lon/Lat
                    all_dF[ii, pp, qq, :] = dF
                    all_dJacobianGWppqq[ii, pp, qq] = dJacobianGWppqq
    return (dFaceArea, dFunIntegral)


def computeAreaIntegralWithGQ(clm, nodes, GN, GW, avg, farea):
    # avg = Boolean flag to take average of the function
    # farea = Boolean flag to compute only the area integral (ignore field)

    # Initialize the area
    dFaceArea = 0.0

    # Initialize the integral
    dFunIntegral = 0.0

    # Set the number of subtriangles
    NST = np.size(nodes, axis=1) - 2

    # Loop over the subtriangles and add up the areas
    # GN, GW = getGaussNodesWeights(order)
    NP = len(GW)

    all_dF = np.zeros(shape=(1, 1, 1, 3), dtype='d')
    all_dJacobianGWppqq = np.zeros(shape=(1, 1, 1), dtype='d')
    if not farea:
        all_dF = np.zeros(shape=(NST, NP, NP, 3), dtype='d')
        all_dJacobianGWppqq = np.zeros(shape=(NST, NP, NP), dtype='d')

    dFaceArea, dFunIntegral = computeAreaIntegralInternal(
        NST, NP, GW, GN, nodes, farea, all_dF, all_dJacobianGWppqq)

    if not farea:
        dFunIntegral = 0.0
        for ii in range(NST):
            # Loop over the quadrature points
            for pp in range(NP):
                for qq in range(NP):
                    dFLonLatRad = trans.computePointCart2LL(
                        all_dF[ii, pp, qq, :])
                    if callable(
                            clm):  # if this is a closed form functional, evaluate directly
                        # print(ii, dF, dFLonLatRad[0], dFLonLatRad[1])
                        thisVar = clm(lon=dFLonLatRad[0], lat=dFLonLatRad[1])
                    else:
                        # Convert to degrees for the SH expansion
                        dFLonLat = 180.0 / mt.pi * dFLonLatRad
                        thisVar = clm.expand(lon=dFLonLat[0], lat=dFLonLat[1])

                    # Sum up the integral of the field
                    dFunIntegral += thisVar * all_dJacobianGWppqq[ii, pp, qq]

    # Compute the cell average
    dFunAverage = dFunIntegral / dFaceArea

    # Return both the integral and the cell average since we computed both
    # Note that when avg = True, dFunAverage = 1.0.
    return [dFunAverage, dFaceArea]

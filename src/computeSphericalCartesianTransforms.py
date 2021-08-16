#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 09:57:18 2019

@author: Jorge Guerra, Vijay Mehadevan
"""

import math as mt
import numpy as np


def computePointCart2LL(pointCoord):

    if len(pointCoord) != 3:
        print('NOT A VALID COORDINATE TRIPLE!')
        return [0.0, 0.0]

    # extract (lon, lat)
    ccoords = pointCoord / np.linalg.norm(pointCoord, ord=2)
    lat = mt.asin(ccoords[2])
    lon = mt.atan2(ccoords[1], ccoords[0])
    if (lon < 0.0):
        lon += 2 * mt.pi
    varLonLat = np.array([lon, lat])

    # OUTPUT IS IN RADIANS
    return varLonLat


def computeCart2LL(cellCoord):
    # Loop over each cell centroid, extract (lon, lat)
    NC = np.size(cellCoord, axis=0)
    varLonLat = np.zeros((NC, 2))
    # print(cellCoord)
    for ii in range(NC):
        ccoords = cellCoord[ii, :] / np.linalg.norm(cellCoord[ii, :], ord=2)
        lat = mt.asin(ccoords[2])
        lon = mt.atan2(ccoords[1], ccoords[0])
        if (lon < 0.0):
            lon += 2 * mt.pi
        varLonLat[ii, :] = [lon, lat]

    # OUTPUT IS IN RADIANS
    return varLonLat


def computeLL2Cart(cellCoord):
    # Loop over the Lon/Lat coordinate array, extract Cartesian coords
    # Input array is [lon, lat, radius]
    NC = np.size(cellCoord, axis=0)
    varCart = np.zeros((NC, 3))
    for ii in range(NC):
        RO = cellCoord[ii, 2]
        lon = cellCoord[ii, 0]
        lat = cellCoord[ii, 1]
        if (lat > 0.5 * mt.pi):
            lat = 0.5 * mt.pi
        if (lat < -0.5 * mt.pi):
            lat = -0.5 * mt.pi
        X = RO * mt.cos(lon) * mt.cos(lat)
        Y = RO * mt.sin(lon) * mt.cos(lat)
        Z = RO * mt.sin(lat)
        RC = mt.sqrt(X**2 + Y**2 + Z**2)
        varCart[ii, :] = [X, Y, Z] / RC

    # INPUT IS IN RADIANS
    return varCart


def computeCentroids(varCon, varCoord):
    # Loop over cells and get centroid vectors
    NC = np.size(varCon, axis=0)
    NP = np.size(varCon, axis=1)
    cellCoord = np.zeros((NC, 3))
    for ii in range(NC):
        # Centroid by averaging corner sphere-center vectors
        centroid = np.array([0.0, 0.0, 0.0])
        for pp in range(NP):
            ndex = varCon[ii, pp] - 1
            centroid += varCoord[:, ndex]

        # Store this centroid
        cellCoord[ii, :] = centroid[:] / NP

    return cellCoord


def computeCentroidsLL(conLon, conLat):
    # Loop over rows of the corner array and get centroid
    NC = np.size(conLon, axis=0)
    NP = np.size(conLon, axis=1)
    cellCoord = np.zeros((NC, 2))
    for ii in range(NC):
        # Centroid by averaging corner sphere-center vectors
        centroid = np.mat([0.0, 0.0])
        for pp in range(NP):
            centroid += [conLon[ii, pp], conLat[ii, pp]]

        centroid *= 1.0 / NP

        # Store this centroid
        cellCoord[ii, :] = centroid

    return cellCoord

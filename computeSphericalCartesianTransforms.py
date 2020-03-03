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
       RO = np.linalg.norm(pointCoord)
       psi = mt.asin(1.0 / RO * pointCoord[2])
       lam = mt.atan2(-pointCoord[0], -pointCoord[1]) + mt.pi
       varLonLat = [lam, psi]
       
       # OUTPUT IS IN RADIANS       
       return np.array(varLonLat)

def computeCart2LL(cellCoord):
       # Loop over each cell centroid, extract (lon, lat)
       NC = np.size(cellCoord, axis=0)
       varLonLat = np.zeros((NC, 2))
       print(cellCoord)
       for ii in range(NC):
              RO = np.linalg.norm(cellCoord[ii,:])
              psi = mt.asin(1.0 / RO * cellCoord[ii,2])
              lam = mt.atan2(-cellCoord[ii,0], -cellCoord[ii,1]) + mt.pi
              varLonLat[ii,:] = [lam, psi]
       
       # OUTPUT IS IN RADIANS       
       return varLonLat

def computeLL2Cart(cellCoord):
       # Loop over the Lon/Lat coordinate array, extract Cartesian coords
       # Input array is [lon, lat, radius]
       NC = np.size(cellCoord, axis=0)
       varCart = np.zeros((NC, 3))
       for ii in range(NC):
              RO = cellCoord[ii,2]
              lon = cellCoord[ii,0]
              lat = cellCoord[ii,1]
              X = RO * mt.cos(lat) * mt.sin(lon)
              Y = RO * mt.cos(lat) * mt.cos(lon)
              Z = RO * mt.sin(lat)
              RC = mt.sqrt(X**2 + Y**2 + Z**2)
              varCart[ii,:] = [X, Y, Z]
              varCart[ii,:] *= 1.0 / RC
       
       # INPUT IS IN RADIANS
       return varCart

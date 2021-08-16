#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:45:49 2019

@author: jeguerra
"""
import numpy as np
import math as mt
from computeSEIntegral import computeSEIntegral
from computeAreaIntegral import computeAreaIntegral

NEL = len(varConGLL)
areaSE = np.zeros((NEL, 1))
areaFV = np.zeros((NEL, 1))

AREASE = 0.0
AREAFV = 0.0
for ii in range(NEL):
    # Compute areas by using the SE method
    cdex = varConGLL[ii, :] - 1
    thisCell = varCoordGLL[:, cdex.astype(int)]
    areaSE[ii], Jacs = computeSEIntegral(thisCell, 4)
    AREASE += areaSE[ii]

    # Compute areas by triangular quadrature on FV mesh
    cdex = varCon[ii, :] - 1
    thisCell = varCoord[:, cdex.astype(int)]
    areaFV[ii] = computeAreaIntegral(None, thisCell, 4, False, True)
    AREAFV += areaFV[ii]

print(4.0 * mt.pi, AREASE, AREAFV)

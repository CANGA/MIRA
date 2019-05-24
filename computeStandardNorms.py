#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 14:28:23 2018

Computes the standard error norms for the regridded and reference target data

@author: jeguerra
"""
import math as mt
import numpy as np
from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

def computeStandardNorms(varConT, varS2T, varST, areaT, jacobiansT, SpectralElement):
       
       # Compute the difference
       varDiff = np.subtract(varS2T, varST)
       varDiff2 = np.power(varDiff, 2)
       varST2 = np.power(varST, 2)
       
       # Compute the necessary integrals
       NT = len(varST)
       L1Den = computeGlobalWeightedIntegral(NT, varConT, varST, areaT, jacobiansT, SpectralElement)
       L2Den = computeGlobalWeightedIntegral(NT, varConT, varST2, areaT, jacobiansT, SpectralElement)
       LinfDen = np.amax(varST)
       
       L1Num = computeGlobalWeightedIntegral(NT, varConT, varDiff, areaT, jacobiansT, SpectralElement)
       L2Num = computeGlobalWeightedIntegral(NT, varConT, varDiff2, areaT, jacobiansT, SpectralElement)
       LinfNum = np.amax(abs(varDiff))
       
       L_1 = L1Num / L1Den
       L_2 = mt.sqrt(L2Num / L2Den)
       L_inf = LinfNum / LinfDen
       
       return L_1, L_2, L_inf
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 09:03:20 2018

Computes the global conservation metric - equation (2) from the document
Returns a single scalar

@author: jeguerra
"""
import numpy as np
from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

def computeGlobalConservation(varConT, varS2T, varST, areaT, jacobiansT, SpectralElement):
       
       # Get the total number of cells/elements
       NT = varConT.shape[0]
       
       diffVarT = np.subtract(varS2T, varST)
       
       # Compute integral of the absolute value of the difference
       L_S2T = computeGlobalWeightedIntegral(NT, varConT, np.abs(diffVarT), areaT, jacobiansT, SpectralElement)
       
       # Compute integral of the (absolute value) target sampled ST data
       L_ST = computeGlobalWeightedIntegral(NT, varConT, np.abs(varST), areaT, jacobiansT, SpectralElement)
       
       # Compute the global conservation metric
       L_g = L_S2T / L_ST
       
       return L_S2T, L_ST, L_g
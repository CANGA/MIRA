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

def computeGlobalConservation(varConS, varConT, varSS, varS2T, varST, areaS, areaT, jacobiansS, jacobiansT, sourceSE, targetSE):
       
       # Get the total number of cells/elements
       NS = varConS.shape[0]
       NT = varConT.shape[0]
       
       # Compute integral of the (absolute value) target sampled ST data
       L_SS = computeGlobalWeightedIntegral(NS, varConS, varSS, areaS, jacobiansS, sourceSE)

       # Compute integral of the absolute value of the difference
       L_S2T = computeGlobalWeightedIntegral(NT, varConT, varS2T, areaT, jacobiansT, targetSE)
       
       # Compute integral of the (absolute value) target sampled ST data
       L_ST = computeGlobalWeightedIntegral(NT, varConT, np.abs(varST), areaT, jacobiansT, targetSE)
       
       # Compute the global conservation metric
       L_g = (L_S2T - L_SS) / L_ST
       # L_g = (L_S2T - L_ST) / L_SS
       
       return L_S2T, L_ST, L_g
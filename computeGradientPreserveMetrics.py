#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 15:04:12 2018

Computes the gradient preservation metrics on the target mesh. Input fields are
lists of numpy arrays containing: regridded and reference gradients and variables
on the target mesh.

@author: TempestGuerra
"""

import numpy as np
import math as mt
from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

def computeGradientPreserveMetrics(varConT, gradsOnTM, varsOnTM, areaT, jacobiansT, SpectralElement):

    # Initialize
    H1 = 0.0
    H1_2 = 0.0
    NT = varConT.shape[0]
         
    # Compute the local errors in the field and gradient
    eK = np.absolute(np.subtract(varsOnTM[1], varsOnTM[0]))
    eGradK = np.linalg.norm(np.subtract(gradsOnTM[1], gradsOnTM[0]), axis=0)
    eK2 = np.multiply(eK, eK)
    eGradK2 = np.multiply(eGradK, eGradK)
    H1num = np.add(eK2, eGradK2)
    
    # Compute some integrals
    varST2 = np.multiply(varsOnTM[0], varsOnTM[0])
    denom = computeGlobalWeightedIntegral(NT, varConT, varST2, areaT, jacobiansT, SpectralElement)
    
    H1num = computeGlobalWeightedIntegral(NT, varConT, H1num, areaT, jacobiansT, SpectralElement)
    H1_2num = computeGlobalWeightedIntegral(NT, varConT, eGradK2, areaT, jacobiansT, SpectralElement)
    
    H1 = mt.sqrt(H1num / denom)
    H1_2 = mt.sqrt(H1_2num / denom)

    return H1, H1_2
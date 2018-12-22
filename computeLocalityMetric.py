#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 10:26:42 2018

Computes the locality metric - equation (5)
Returns an array of length 'number of elements' corresponding to locality
measure at each target element (grid cell)

@author: jeguerra
"""

import numpy as np
from computeAreaWeight import computeAreaWeight
from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

def computeLocalityMetric(varSS, varST, varConS, varCoordS, varConT, varCoordT):
       
       # Get the length of the arrays
       NS = len(varSS)
       NT = len(varST)
       
       # Precompute the area weights and then look them up in the integrals below
       areas = np.zeros((NT,1))
       for ii in range(NT):
              areas[ii] = computeAreaWeight(varCoordT, varConT[ii, :])
       
       # Initialize the columns of the R*E operator (see eq. 4 for E as a scaling operator)
       varDiffT = np.zeros(NT,1)
       NumSum = 0.0
       DemSum = 0.0
       
       for ii in range(NS):
              for jj in range(NT):
                     # Compute distance from this target location to the source location
                     dij = 0.0
                     # Get the regridded "delta field" E_i^s (NOT POSSIBLE WITHOUT R MAP)
                     Esi = 0.0
                     # Compute the difference over the target mesh
                     varDiffT[ii] = dij * abs(Esi - varST)
                     
              # Compute integral of the target locality data
              L_DF = computeGlobalWeightedIntegral(NT, varConT, varCoordT, varDiffT, areas)
                     
              # Compute integral of the target sampled ST data
              L_ST = computeGlobalWeightedIntegral(NT, varConT, varCoordT, varST, areas)
              
              NumSum += L_DF
              DemSum += L_ST
              
       L_l = NumSum / DemSum
        
       return L_l
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 09:03:20 2018

Computes the global conservation metric - equation (2) from the document
Returns a single scalar

@author: jeguerra
"""

from computeGlobalWeightedIntegral import computeGlobalWeightedIntegral

def computeGlobalConservation(varSS, varS2T, varST, areaS, areaT):
       
       NS = len(varSS)
       NT = len(varST)
       
       # Compute integral of the source sampled SS data
       L_SS = computeGlobalWeightedIntegral(NS, varSS, areaS)
       
       # Compute integral of the source to target remapped S2T data
       L_S2T = computeGlobalWeightedIntegral(NT, varS2T, areaT)
       
       # Compute integral of the target sampled ST data
       L_ST = computeGlobalWeightedIntegral(NT, varST, areaT)
       
       # Compute the global conservation metric
       L_g = (L_S2T - L_SS) / L_ST
       
       return L_g[0]
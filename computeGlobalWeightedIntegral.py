#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 12:24:42 2018

Computes a global area integral based on Exodus and field variable data. The
various data arrays MUST have matching sizes. This is checked here.

@author: jeguerra
"""
import numpy as np

def computeGlobalWeightedIntegral(NEL, varCon, varF, areas, jacobians, SpectralElement):
       # NEL = total number of elements/cells
       
       # Loop over each element and compute the sum
       INT = 0.0
       if not SpectralElement:
              INT = varF.dot(areas)

       else:
             for ii in range(NEL):
                    gdex = varCon[ii,:] - 1
                    INT += np.dot(jacobians[ii,:], varF[gdex.astype(int)])
              
       return INT

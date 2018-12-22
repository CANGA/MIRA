#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 12:24:42 2018

Computes a global area integral based on Exodus and field variable data. The
various data arrays MUST have matching sizes. This is checked here.

@author: jeguerra
"""

def computeGlobalWeightedIntegral(NEL, varF, areas):
       
       try:
              # Check the number of elements (grid cells) against field variable
              if len(varF) != NEL:
                     print('Field variable array length does not match number of elements in "computeGlobalWeightedIntegral"')
       except ValueError:
              print('CHECK DATA ARRAY SIZES FOR INTEGRAL COMPUTATIONS.')
              
       # Loop over each element and compute the sum
       INT = 0.0
       for ii in range(NEL):
              INT += areas[ii] * varF[ii]
              
       return INT
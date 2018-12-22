#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 08:15:32 2018

Computes the global extrema metrics for regridded and reference target data

@author: jeguerra
"""

import numpy as np

def computeGlobalExtremaMetrics(varS2T, varST):
       
       L_min = (np.amin(varST) - np.amin(varS2T)) / np.min(abs(varST))
       L_max = (np.amax(varS2T) - np.amax(varST)) / np.max(abs(varST))
       
       return L_min, L_max
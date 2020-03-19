#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 19 08:15:32 2018

Computes the global extrema metrics for regridded and reference target data

@author: jeguerra
"""

import numpy as np

def computeGlobalExtremaMetrics(varS2T, varST):
       
       Lden = (np.max(abs(varST)) - np.min(abs(varST)))
       L_min = np.abs(np.amin(varST) - np.amin(varS2T)) / Lden
       L_max = np.abs(np.amax(varS2T) - np.amax(varST)) / Lden
       
       return L_min, L_max
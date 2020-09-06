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
    L_min = min(np.amin(varS2T) - np.amin(varST), 0.0) / Lden  # < 0 indicates failure
    L_max = max(np.amax(varS2T) - np.amax(varST), 0.0) / Lden  # > 0 indicates failure

    return L_min, L_max

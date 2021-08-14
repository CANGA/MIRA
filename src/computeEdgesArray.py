#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 09:27:36 2019

 Computes the local NP X 2 array of edges based on nodal connectivity

@author: jeguerra
"""

import numpy as np


def computeEdgesArray(NP, connect):

    edex = np.zeros(2 * NP, dtype=int)
    jj = 1
    for ii in range(NP - 1):
        edex[jj] = ii + 1
        edex[jj + 1] = ii + 1
        jj += 2

    edges = np.reshape(connect[edex], (NP, 2))

    return edges

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 25 16:28:19 2019

@author: TempestGuerra
"""

import numpy as np


def computeCentroid(NP, cell):
    # Centroid by averaging corner sphere-center vectors
    centroid = np.mat([0.0, 0.0, 0.0])
    for pp in range(NP):
        centroid += cell[:, pp]

    centroid *= 1.0 / NP

    # Renormalize the centroid vector
    RO = np.linalg.norm(centroid)
    centroid *= 1.0 / RO

    return centroid

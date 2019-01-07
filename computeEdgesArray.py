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
              
       """
       if NP == 4:
              edex = [0, 1, 1, 2, 2, 3, 3, 0]
              
              edges = np.array([[connect.flat[edex[0]], connect.flat[edex[1]]], \
                                [connect.flat[edex[2]], connect.flat[edex[3]]], \
                                [connect.flat[edex[4]], connect.flat[edex[5]]], \
                                [connect.flat[edex[6]], connect.flat[edex[7]]]])
       elif NP == 3:
              edex = [0, 1, 1, 2, 2, 0]
                     
              edges = np.array([[connect.flat[edex[0]], connect.flat[edex[1]]], \
                                [connect.flat[edex[2]], connect.flat[edex[3]]], \
                                [connect.flat[edex[4]], connect.flat[edex[5]]]])
       else:
              print('ONLY TRIANGLES OR QUADRILATERALS SUPPORTED!')
       """
       
       return edges

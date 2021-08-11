#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(
    name='CANGA-ROO',
    version='1.0.0',
    packages=find_packages(include=['src']),
    install_requires=[
        'pyshtools>=4.8',
        'scipy>=1.7',
        'plotly>=5.1.0',
        'pandas>=1.3.0',
        'matplotlib>=3.4.0',
        'netCDF4>=1.5.7',
        'numba>=0.53.1',
        'rasterio>=1.2.6'
    ]
)

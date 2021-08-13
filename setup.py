#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

setup(
    name='CANGA-ROO',
    version='1.0.0',
    packages=find_packages(include=['src']),
    install_requires=[
        'pyshtools>=4.2',
        'scipy>=1.0.1',
        'plotly>=4.0.0',
        'pandas>=0.23.0',
        'matplotlib>=2.2.0',
        'netCDF4>=1.4.0',
        'rasterio>=1.0.0',
        'numba>=0.48.0',
    ]
)

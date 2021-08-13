#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='CANGA-ROO',
    version='1.0.0',
    description="Coupling Approaches for Next-Generation Architectures (CANGA) -- Remapping Offline-Online (ROO) Intercomparison Package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/CANGA/Remapping-Intercomparison",
    project_urls={
        "Homepage": "https://www.canga-scidac.org/"
        "Bug Tracker": "https://github.com/CANGA/Remapping-Intercomparison/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD 3-clause License",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(include=['src']),
    install_requires=[
        'pyshtools>=4.2',
        'scipy>=1.0.1',
        'plotly>=4.0.0',
        'pandas>=0.23.0',
        'matplotlib>=2.2.0',
        'netCDF4>=1.4.0',
        'rasterio>=1.0.0',
        'numba>=0.48.0'
    ],
    python_requires=">=3.6"
)

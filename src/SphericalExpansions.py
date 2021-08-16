#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 09:25:34 2019

Code for preprocessing spherical harmonic expansion of global topography,
total precipitable water composite, and cloud fraction. Topographic and TPW data
are at 0.25 degree resolution.

@author: jeguerra

"""
# %%
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.optimize import curve_fit
from scipy.optimize import newton
from mpl_toolkits.mplot3d import Axes3D
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
import pyshtools
import rasterio as rio

# %% PROCESS THE TPW DATA

# Make a list of input files (4 samples)
dataDir = 'TPW_DATA/'
tpw_files = ['comp20190209.160000.nc', 'comp20190210.180000.nc',
             'comp20190211.080000.nc', 'comp20190211.210000.nc']

tpw_coeffs = []

# Loop over the input files
NF = len(tpw_files)
for ii in range(NF):
    # Read in the tpw data
    tpw_fid = Dataset(dataDir + tpw_files[ii])
    lon = tpw_fid.variables['lonArr'][:]
    lat = tpw_fid.variables['latArr'][:]
    tpw = tpw_fid.variables['tpwGrid'][:]

    # Fix TPW data so that lon data is 2X lat data
    tpw = tpw[1:, :]

    # Perform the SH expansion
    coeffs = pyshtools.expand.SHExpandDH(tpw, sampling=2)
    tpw_coeffs.append(coeffs)

# %% Average the 4 coefficients
avgCoeffsTPW = np.average(tpw_coeffs, axis=0)

# %% Compute the spectrum of TPW data
powerTPW = pyshtools.spectralanalysis.spectrum(avgCoeffsTPW, unit='per_l')
degsTPW = np.arange(avgCoeffsTPW.shape[1])

# %% PROCESS THE CLOUD FRACTION DATA

# Make a list of input files (4 samples in grayscale)
dataDir = 'CloudFraction_DATA/'
cf_files = ['MODAL2_E_CLD_FR_2018-11-17_gs_3600x1800.TIFF',
            'MODAL2_E_CLD_FR_2018-07-20_gs_3600x1800.TIFF',
            'MODAL2_E_CLD_FR_2018-09-22_gs_3600x1800.TIFF',
            'MODAL2_E_CLD_FR_2018-11-17_gs_3600x1800.TIFF']

cf_coeffs = []

# Loop over the input files
NF = len(cf_files)
for ii in range(NF):
    # Read in the cloud fraction data (raster data with max = 255)
    cfrac = rio.open(dataDir + cf_files[ii])
    cf_data = 1.0 / 255.0 * cfrac.read(1)

    # Perform the SH expansion
    coeffs = pyshtools.expand.SHExpandDH(cf_data, sampling=2)
    cf_coeffs.append(coeffs)

# %% Average the 4 coefficients
avgCoeffsCF = np.average(cf_coeffs, axis=0)

# %% Compute the spectrum of Cloud Fraction data
powerCF = pyshtools.spectralanalysis.spectrum(avgCoeffsCF, unit='per_l')
degsCF = np.arange(avgCoeffsCF.shape[1])

# %% PROCESS THE TOPOGRAPHY DATA
# Test data from the shtools example
#clm, lmax = pyshtools.shio.shread('srtmp300.msl')
#topo = pyshtools.expand.MakeGridDH(clm, sampling=2)

# Read in the NASA data from https://www.ngdc.noaa.gov/mgg/global/global.html
dataDir = 'Topography_DATA/'
topo_file = 'ETOPO1_Ice_g_gmt4.grd'
topo_fid = Dataset(dataDir + topo_file)
topo = topo_fid.variables['z'][:]

# Fix topo data for the SH expansion
topo = topo[1:, 1:]

# Get the coefficients
coeffsTOPO = pyshtools.expand.SHExpandDH(topo, sampling=2)

# %% Compute the spectrum of topography data
powerTOPO = pyshtools.spectralanalysis.spectrum(coeffsTOPO, unit='per_l')
degsTOPO = np.arange(coeffsTOPO.shape[1])

# %% Save the coefficients
saveOut = {
    'avgCoeffsTPW': avgCoeffsTPW,
    'avgCoeffsCF': avgCoeffsCF,
    'coeffsTOPO': coeffsTOPO}
np.savez('SHCoefficientsRaw_TPW_CFR_TOPO', **saveOut)

# %% Do a curve fit to find the monomial relation


def func(x, a, b, c):
    return a * np.power(x, b) + c


# Convert to amplitude from power spectra
#ampTPW = np.power(powerTPW, 0.5)
#ampCF = np.power(powerCF, 0.5)
#ampTOPO = np.power(powerTOPO, 0.5)
ampTPW = powerTPW
ampCF = powerCF
ampTOPO = powerTOPO
endPoint = 1001

bounds1 = ([0.0, -50.0, -np.inf], [np.inf, 50.0, np.inf])
bounds2 = ([0.0, -10.0, -10.0], [np.inf, 0.0, 10.0])
# Fit the TPW data (INDEXED TO EXTRACT SUBSET THAT IS BETTER FOR FIT)
poptTPW1, pcovTPW1 = curve_fit(func, degsTPW[1:5], ampTPW[1:5], bounds=bounds1)
poptTPW2, pcovTPW2 = curve_fit(
    func, degsTPW[20:endPoint], ampTPW[20:endPoint], bounds=bounds2)
# Fit the Cloud Fraction data (INDEXED TO EXTRACT SUBSET THAT IS BETTER
# FOR FIT)
poptCF1, pcovCF1 = curve_fit(func, degsCF[1:21], ampCF[1:21], bounds=bounds1)
poptCF2, pcovCF2 = curve_fit(
    func, degsCF[20:endPoint], ampCF[20:endPoint], bounds=bounds2)
# Fit the topography data (INDEXED TO EXTRACT SUBSET THAT IS BETTER FOR FIT)
bounds1 = ([0.0, -100.0, -np.inf], [np.inf, 50.0, np.inf])
bounds2 = ([0.0, -np.inf, -np.inf], [np.inf, 0.0, np.inf])
poptTOPO1, pcovTOPO1 = curve_fit(
    func, degsTOPO[1:5], ampTOPO[1:5], bounds=bounds1)
poptTOPO2, pcovTOPO2 = curve_fit(
    func, degsTOPO[20:endPoint], ampTOPO[20:endPoint], bounds=bounds2)

# %% Find the intersection point for each spectrum pair


def intFunc(x, a1, b1, c1, a2, b2, c2):
    return a1 * np.power(x, b1) - a2 * np.power(x, b2) + (c1 - c2)


params = np.append(poptTPW1, poptTPW2)
degIntxTPW = newton(intFunc, x0=2.0, args=(params))

params = np.append(poptCF1, poptCF2)
degIntxCF = newton(intFunc, x0=5.0, args=(params))

params = np.append(poptTOPO1, poptTOPO2)
degIntxTOPO = newton(intFunc, x0=4.0, args=(params))

print(degIntxTPW, degIntxCF, degIntxTOPO)

# %% PLOT THE TPW AND TOPO DATA SPECTRA
with plt.style.context(('ggplot')):
    fig, (ax0, ax1, ax2) = plt.subplots(
        nrows=3, figsize=(12, 10), tight_layout=True)

    # Plot the TPW spectrum
    fitTPW1 = func(degsTPW, *poptTPW1)
    fitTPW2 = func(degsTPW, *poptTPW2)
    ax0.plot(degsTPW, ampTPW, 'k')
    ax0.plot(degsTPW, fitTPW1, 'b--',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(poptTPW1))
    ax0.plot(degsTPW, fitTPW2, 'g--',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(poptTPW2))
    ax0.plot([degIntxTPW], [func(degIntxTPW, *poptTPW1)],
             marker='D', markersize=5, label='intersect = %5.3f' % degIntxTPW)
    ax0.set_title('Total Precipitable Water - SH Amplitude Spectrum')
    ax0.set(yscale='log', xscale='log', ylabel='Avg. Coefficient Amplitude')
    ax0.legend()
    ax0.grid(b=True, which='both', axis='both')
    # Plot the Cloud Fraction spectrum
    fitCF1 = func(degsCF, *poptCF1)
    fitCF2 = func(degsCF, *poptCF2)
    ax1.plot(degsCF, ampCF, 'k')
    ax1.plot(degsCF, fitCF1, 'b--',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(poptCF1))
    ax1.plot(degsCF, fitCF2, 'g--',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(poptCF2))
    ax1.plot([degIntxCF], [func(degIntxCF, *poptCF1)],
             marker='D', markersize=5, label='intersect = %5.3f' % degIntxCF)
    ax1.set_title('Global Cloud Fraction - SH Amplitude Spectrum')
    ax1.set(yscale='log', xscale='log', ylabel='Avg. Coefficient Amplitude')
    ax1.legend()
    ax1.grid(b=True, which='both', axis='both')
    # Plot the Topography spectrum
    fitTOPO1 = func(degsTOPO, *poptTOPO1)
    fitTOPO2 = func(degsTOPO, *poptTOPO2)
    ax2.plot(degsTOPO, ampTOPO, 'k')
    ax2.plot(degsTOPO, fitTOPO1, 'b--',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(poptTOPO1))
    ax2.plot(degsTOPO, fitTOPO2, 'g--',
             label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(poptTOPO2))
    ax2.plot([degIntxTOPO], [func(degIntxTOPO, *poptTOPO1)],
             marker='D', markersize=5, label='intersect = %5.3f' % degIntxTOPO)
    ax2.set_title('Global Topography Data - SH Amplitude Spectrum')
    ax2.set(
        yscale='log',
        xscale='log',
        xlabel='Spherical harmonic degree',
        ylabel='Avg. Coefficient Amplitude')
    ax2.legend()
    ax2.grid(b=True, which='both', axis='both')
plt.show()

plt.savefig("TPW_Topo_CloudFrac_SHpowerSpectra.png")

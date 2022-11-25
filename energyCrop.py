# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 20:17:49 2022

@author: vicke
"""

import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import netCDF4 as nc
from netCDF4 import Dataset


def globalarea(reso):
    lon=int(360/reso)
    lat=int(180/reso)
    lat_area = np.zeros(lat);RE = 6371.393 #earth redius 6371.393 k
    for i in np.arange(lat):
      lat_area[i] = abs(2*np.pi*(RE**2)*(np.cos(np.pi*(i+1)/lat)-np.cos(np.pi*i/lat)))/lon
    area = np.zeros((lat,lon))
    for i in np.arange(lon):
      area[:,i] = lat_area
    return(area)

area = globalarea(0.01)
chinaArea = area[10500: 14500, 7000: 14000]

# ------------ 1. energy crop dataset ----------------
# read ecrop nc file
ecrop = Dataset('../inputData/map/eCrop/Bioenergy_crop_yields.nc','r')
bestCrop = np.array(ecrop.variables['Best_crop'])[70:150,500:640] # 0.5 degree

# reclass dataset
bestCrop_1km = np.zeros([4000,7000])
for i in np.arange(bestCrop.shape[0]):
    for j in np.arange(bestCrop.shape[1]):
        bestCrop_1km[i*50:(i+1)*50,j*50:(j+1)*50]=bestCrop[i,j]

bestCrop_1km = bestCrop_1km * 100 * chinaArea  # 1t/hm = 100 tonne/km
bestCrop_1km = np.nan_to_num(bestCrop_1km, nan=0)

# ------------ 2. energy crop on marginalLand ---------
ncfile = nc.Dataset('../outputData/bioenergy.nc','a',format='NETCDF4')
land1 = np.array(ncfile.variables['land_margin'])
land2 = np.array(ncfile.variables['land_abandon'])

eCropLand1 = land1 * bestCrop_1km # max
eCropLand2 = land2 * bestCrop_1km # min

# write results in nc file

# create variable
nc_eCrop1 = ncfile.createVariable('eCrop_abandon', np.float32, ('lat','lon'))
nc_eCrop1.units = 'tonne/km2'
nc_eCrop1.standard_name = 'energy_crop_on_marginal_land_min'
nc_eCrop1[:,:] = eCropLand2

# create variable
nc_eCrop2 = ncfile.createVariable('eCrop_margin', np.float32, ('lat','lon'))
nc_eCrop2.units = 'tonne/km2'
nc_eCrop2.standard_name = 'energy_crop_on_marginal_land_max'
nc_eCrop2[:,:] = eCropLand1

# close file
ncfile.close()

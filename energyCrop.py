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
miscanthus = np.array(ecrop.variables['Miscanthus'])[70:150,500:640]
eucalypt = np.array(ecrop.variables['Eucalypt'])[70:150,500:640]
poplar = np.array(ecrop.variables['Poplar'])[70:150,500:640]
switchgrass = np.array(ecrop.variables['Switchgrass'])[70:150,500:640]
willow = np.array(ecrop.variables['Willow'])[70:150,500:640]

energyCropList = [miscanthus, eucalypt, poplar, switchgrass, willow]

gridNum = 4000 * 7000
eCropsArray= np.zeros([gridNum, 8])
k = 0

energyCropList_1km = []
for eArray in energyCropList:
    eArray_1km = np.zeros([4000, 7000])
    for i in np.arange(eArray.shape[0]):
        for j in np.arange(eArray.shape[1]):
            eArray_1km[i*50:(i+1)*50, j*50:(j+1)*50] = eArray[i,j]
    eArray_1km = eArray_1km * 100 * chinaArea
    eArray_1km = np.nan_to_num(eArray_1km, nan=0)
    energyCropList_1km.append(eArray_1km)  
    eCropsArray[:,k] = eArray_1km.reshape(gridNum)
    k = k + 1

print(k)    
eCropsArray[:,5] = np.nanmax(eCropsArray[:,[0,1,2,3,4]],axis=1)
eCrops_best = eCropsArray[:,5].reshape(4000,7000)
energyCropList_1km.append(eCrops_best)

# ------------ 2. energy crop on marginalLand ---------
ncfile = nc.Dataset('../outputData/NetCDF/BioRes_breakdown.nc','a',format='NETCDF4')
land1 = np.array(ncfile.variables['land_margin'])
land2 = np.array(ncfile.variables['land_abandon'])

eCropName = ['miscanthus', 'eucalypt', 'poplar', 'switchgrass', 'willow', 'best']
name_i = 0
for eArray_1km in energyCropList_1km[:]:
    
    name = 'eCrops_' + eCropName[name_i]
    
    eCrop_margin = eArray_1km * land1
    eCrop_abandon = eArray_1km * land2
    
    # create variable
    nc_eMargin = ncfile.createVariable(name + "_margin", np.float32, ('lat','lon'))
    nc_eMargin.units = 'tonne'
    nc_eMargin.standard_name = name + "_margin"
    nc_eMargin[:,:] = eCrop_margin
    
    # create variable
    nc_eAbandon = ncfile.createVariable(name + "_abandon", np.float32, ('lat','lon'))
    nc_eAbandon.units = 'tonne'
    nc_eAbandon.standard_name = name + "_abandon"
    nc_eAbandon[:,:] = eCrop_abandon
    
    name_i = name_i+1

ncfile.close()
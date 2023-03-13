# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 19:08:29 2022

@author: vicke
"""

import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import netCDF4 as nc

# -------------- 1. marginal land by type ------------------
# (1) read input dataset
ds = gdal.Open('../inputData/map/LandUseType/landUse2018.tif')
landUse2018 = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/excludeArea/protectedArea2.tif')
protectLand = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/excludeArea/Slope.tif')
slope = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
chinaProv = ds.GetRasterBand(1).ReadAsArray()

# (2) filter land by type
marginLand = landUse2018.copy()
marginLand[(marginLand==22)|(marginLand==23)|(marginLand==31)|(marginLand==32)|(marginLand==33)|
           (marginLand==45)|(marginLand==46)|(marginLand==63)|(marginLand==65)] = 1
marginLand[marginLand!=1] = 0

# (3) exclude protect area
protectLand[protectLand!=1] = 0 # cannot select protect area
exclude_protect = 1-protectLand

# (4) exclude slope larger than 25
slope[slope<=25] = 1
slope[slope>25] = 0 # cannot select slope larger than 25

# (5) exclude graze area
grazeLand = landUse2018.copy()
grazeLand[(grazeLand==31)|(grazeLand==32)] = 1
grazeLand[grazeLand !=1 ] = 0
prov_grazeLand = chinaProv.copy()
prov_grazeLand[(prov_grazeLand==31)|(prov_grazeLand==22)|(prov_grazeLand==20)|(prov_grazeLand==21)|(prov_grazeLand==32)] = 1
prov_grazeLand[prov_grazeLand != 1] = 0
exclude_grazeLand = grazeLand * prov_grazeLand
exclude_grazeLand = 1 - grazeLand # cannot select grazeLand

# (6) all criteria
marginLand = marginLand * exclude_protect * slope * exclude_grazeLand

# ------------- 2. marginal land by time sequential ------
ds = gdal.Open('../inputData/map/LandUseType/landUse1980.tif')
landUse1980 = ds.GetRasterBand(1).ReadAsArray()
crop_land1980 = landUse1980.copy()
crop_land1980[(crop_land1980 == 11) | (crop_land1980 == 12)] = 1
crop_land1980[crop_land1980 != 1] = 0

crop_land2018 = landUse2018.copy()
crop_land2018[(crop_land2018 == 11) | (crop_land2018 == 12)] = 1
crop_land2018[crop_land2018 != 1] = 0 

abandonCrop = (1-crop_land2018) - (1-crop_land1980)
abandonCrop[abandonCrop <1 ] = 0

abandonCrop = abandonCrop * exclude_protect

# ----------- 3. write as nc file ----------------------
# create dimension; create variable; fill dimension; fill variables
ncfile = nc.Dataset('../outputData/NetCDF/BioRes_breakdown.nc','w',format='NETCDF4')
ncfile.createDimension('lat', 4000) # latitude axis
ncfile.createDimension('lon', 7000)  # longitude axis

# create variable
land1 = ncfile.createVariable('land_abandon', np.float32, ('lat','lon'))
land1.units = 'km2'
land1.standard_name = 'marginal_land_by_time_sequential'
land1[:,:] = abandonCrop

land2 = ncfile.createVariable('land_margin', np.float32, ('lat','lon'))
land2.units = 'km2'
land2.standard_name = 'marginal_land_by_type'
land2[:,:] = marginLand

ncfile.close()
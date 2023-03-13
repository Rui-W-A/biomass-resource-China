# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 21:09:04 2022

@author: vicke
"""

import numpy as np
import netCDF4 as nc
import pandas as pd
from osgeo import gdal

def writeTif_wgs84(pathOut,dataMatrix):
    out_ds=gdal.GetDriverByName('GTiff').Create(
        pathOut,
        dataMatrix.shape[1], # 矩阵的列数col
        dataMatrix.shape[0], # 矩阵的行数row
        1,
        gdal.GDT_Float64,
        options=['COMPRESS=LZW'])
    ds=gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
    out_ds.SetProjection(ds.GetProjection())
    res1 = ds.GetGeoTransform()[1];res2 = ds.GetGeoTransform()[-1];
    GetGeoTransform = [70,res1,0,55,0,res2]
    out_ds.SetGeoTransform(GetGeoTransform)
    out_band = out_ds.GetRasterBand(1).WriteArray(dataMatrix)
    del out_ds,out_band,ds

# read BioRes_breakdown
ncfile = nc.Dataset('../outputData/NetCDF/BioRes_breakdown.nc','r',format='NETCDF4')

# write agricultural residues as Geotiff
agri_columns = ['rice', 'maize', 'wheat', 'other grain', 'cotton', 'canola', 'peanuts', 'soybeans', 'potatoes', 'all']
for i in agri_columns:
    name = 'agri_' + i
    tmp = np.array(ncfile.variables[name])
    writeTif_wgs84('../outputData/Geotiff/' + name +'.tif', tmp)

# write forestry residues as Geotiff
fores_columns = ['woodprocess_merchant', 'woodprocess_fire', 'woodprocess_log', 'woodprocess_bamboo',
             'logging_fuel', 'logging_timber', 'logging_protect', 'logging_economic',
             'logging_specialpurpose', 'logging_shrub', 'logging_sparse', 'all']
for i in fores_columns[:]:
    name = 'fores_'+ i
    tmp = np.array(ncfile.variables[name])
    writeTif_wgs84('../outputData/Geotiff/' + name +'.tif', tmp)

# write energy crops as Geotiff
eCropName = ['best', 'miscanthus', 'eucalypt', 'poplar', 'switchgrass', 'willow']
for i in eCropName:
    margin_i = 'eCrops_' + i + '_margin'
    tmp_margin = np.array(ncfile.variables[margin_i])
    writeTif_wgs84('../outputData/Geotiff/' + margin_i + '.tif', tmp_margin)
    
    abandon_i = 'eCrops_' + i + '_abandon'
    tmp_abandon = np.array(ncfile.variables[abandon_i])
    writeTif_wgs84('../outputData/Geotiff/' + abandon_i + '.tif', tmp_abandon)

# write scenario data as Geotiff
sceName = ['s1_AFE_max', 's2_AFE_margin', 's3_AFE_abandon', 's4_AF_max','s5_AF_min']
nc_aggre = nc.Dataset('../outputData/NetCDF/BioEnergy_scenario.nc', 'r', format='NETCDF4')
for i in sceName:
    tmp = np.array(nc_aggre.variables[i])
    writeTif_wgs84('../outputData/Geotiff/'+ 'bioenergy_' + i + '.tif', tmp)
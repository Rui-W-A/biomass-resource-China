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
    # 左上角x坐标， 水平分辨率，旋转参数， 左上角y坐标，旋转参数，竖直分辨率。
    GetGeoTransform = [70,res1,0,55,0,res2]
    out_ds.SetGeoTransform(GetGeoTransform)
    out_band = out_ds.GetRasterBand(1).WriteArray(dataMatrix)
    del out_ds,out_band,ds

ncfile = nc.Dataset('../outputData/bioenergy.nc','a',format='NETCDF4')
agri = np.array(ncfile.variables['agri_Res'])
fores = np.array(ncfile.variables['fores_Res'])
eCrop_abandon = np.array(ncfile.variables['eCrop_abandon'])
eCrop_margin = np.array(ncfile.variables['eCrop_margin'])
bio_s1 = np.array(ncfile.variables['AFE_max'])
bio_s2 = np.array(ncfile.variables['AFE_min'])
bio_s3 = np.array(ncfile.variables['AF_max'])
bio_s4 = np.array(ncfile.variables['AF_min'])

writeTif_wgs84('../outputData/agri.tif', agri)
writeTif_wgs84('../outputData/fores.tif', fores)
writeTif_wgs84('../outputData/eCrop1_margin.tif', eCrop_margin)
writeTif_wgs84('../outputData/eCrop2_abandon.tif', eCrop_abandon)
writeTif_wgs84('../outputData/bioenergy_s1.tif', bio_s1)
writeTif_wgs84('../outputData/bioenergy_s2.tif', bio_s2)
writeTif_wgs84('../outputData/bioenergy_s3.tif', bio_s3)
writeTif_wgs84('../outputData/bioenergy_s4.tif', bio_s4)

# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:48:36 2022

@author: vicke
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
from osgeo import gdal

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

def zonalByProvin_Table(chinaProvinTif,chinaTable,temData):
    sumTable=chinaTable.copy()
    # sumTable=sumTable.set_index('OBJECTID')
    for i in range(1,35,1):
        tem_energyProd=np.copy(temData)
        # 各省份掩膜
        tem_chinaTif=np.copy(chinaProvinTif)
        tem_chinaTif[tem_chinaTif!=i]=0
        tem_chinaTif[tem_chinaTif!=0]=1
        # 用掩膜乘以产量数据
        tem_provinTif=tem_energyProd*tem_chinaTif
        sumResult=np.nansum(tem_provinTif)
        sumTable.loc[i,'sum']=sumResult
    return sumTable

# read data
ncfile = nc.Dataset('../outputData/bioenergy.nc','r',format='NETCDF4')
land_min = np.array(ncfile.variables['land_abandon']) * chinaArea
land_max = np.array(ncfile.variables['land_margin']) * chinaArea
agri = np.array(ncfile.variables['agri_Res'])
fores = np.array(ncfile.variables['fores_Res'])
eCrop_min = np.array(ncfile.variables['eCrop_abandon']) 
eCrop_max = np.array(ncfile.variables['eCrop_margin'])
bio_af_min = np.array(ncfile.variables['AF_min'])
bio_af_max = np.array(ncfile.variables['AF_max'])
bio_afe_min = np.array(ncfile.variables['AFE_min'])
bio_afe_max = np.array(ncfile.variables['AFE_max'])
ncfile.close()

# read china province
chinaTable = pd.read_excel('../inputData/map/ChinaMask/ChinaTable.xlsx', index_col='OBJECTID')
ds = gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
chinaProv = ds.GetRasterBand(1).ReadAsArray()

# zonal by province
prov_land_min = zonalByProvin_Table(chinaProv, chinaTable, land_min)
prov_land_max = zonalByProvin_Table(chinaProv, chinaTable, land_max)
prov_agri = zonalByProvin_Table(chinaProv, chinaTable, agri)
prov_fores = zonalByProvin_Table(chinaProv, chinaTable, fores)
prov_eCrop_min = zonalByProvin_Table(chinaProv, chinaTable, eCrop_min)
prov_eCrop_max = zonalByProvin_Table(chinaProv, chinaTable, eCrop_max)
prov_bio_af_min = zonalByProvin_Table(chinaProv, chinaTable, bio_af_min)
prov_bio_af_max = zonalByProvin_Table(chinaProv, chinaTable, bio_af_max)
prov_bio_afe_min = zonalByProvin_Table(chinaProv, chinaTable, bio_afe_min)
prov_bio_afe_max = zonalByProvin_Table(chinaProv, chinaTable, bio_afe_max)

prov_out = chinaTable.copy()
prov_out.loc[:,'land_abandon_km2'] = prov_land_min.loc[:,'sum']
prov_out.loc[:,'land_margin_km2'] = prov_land_max.loc[:,'sum']
prov_out.loc[:,'agri_Mt'] = prov_agri.loc[:,'sum'] / pow(10,6)
prov_out.loc[:,'fores_Mt'] = prov_fores.loc[:,'sum'] / pow(10,6)
prov_out.loc[:,'eCrop_abandon_Mt'] = prov_eCrop_min.loc[:,'sum'] / pow(10,6)
prov_out.loc[:,'eCrop_margin_Mt'] = prov_eCrop_max.loc[:,'sum'] / pow(10,6)
prov_out.loc[:,'bio_AF_min_EJ'] = prov_bio_af_min.loc[:,'sum'] / pow(10,9)
prov_out.loc[:,'bio_AF_max_EJ'] = prov_bio_af_max.loc[:,'sum'] / pow(10,9)
prov_out.loc[:,'bio_AFE_min_EJ'] = prov_bio_afe_min.loc[:,'sum'] / pow(10,9)
prov_out.loc[:,'bio_AFE_max_EJ'] = prov_bio_afe_max.loc[:,'sum'] / pow(10,9)

prov_out.to_csv('../outputData/prov_bioenergy.csv')






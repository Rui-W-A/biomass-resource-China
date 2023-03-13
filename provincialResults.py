# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:48:36 2022

@author: vicke
"""

import pandas as pd
import numpy as np
import netCDF4 as nc
from osgeo import gdal

def zonalByProvin_Table(chinaProvinTif,chinaTable,temData):
    sumTable=chinaTable.copy()
    # sumTable=sumTable.set_index('OBJECTID')
    for i in range(1,35,1):
        tem_energyProd=np.copy(temData)
        # mask for each province
        tem_chinaTif=np.copy(chinaProvinTif)
        tem_chinaTif[tem_chinaTif!=i]=0
        tem_chinaTif[tem_chinaTif!=0]=1
        # multiply mask with production
        tem_provinTif=tem_energyProd*tem_chinaTif
        sumResult=np.nansum(tem_provinTif)
        sumTable.loc[i,'sum']=sumResult
    return sumTable

# China data
chinaTable = pd.read_excel('../inputData/map/ChinaMask/ChinaTable.xlsx', index_col='OBJECTID')
prov_out_land = chinaTable.copy()
prov_out_agri = chinaTable.copy()
prov_out_fores = chinaTable.copy()
prov_out_eCrops = chinaTable.copy()
prov_out_bioScen = chinaTable.copy()
ds = gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
chinaProv = ds.GetRasterBand(1).ReadAsArray()

# read netCDF dataset
# BioRes_breakdown
ncfile = nc.Dataset('../outputData/NetCDF/BioRes_breakdown.nc','r',format='NETCDF4')
ncfile_keys = list(ncfile.variables.keys())
for key in ncfile_keys:
    tmp = np.array(ncfile.variables[key])
    prov_tmp = zonalByProvin_Table(chinaProv, chinaTable, tmp)
    if 'land' in key:  # km2
        print(key)
        colName = key + '_km2'
        prov_out_land[colName] = prov_tmp.loc[:,'sum']
        
    elif 'agri' in key:
        print(key)
        colName = key + '_Mt'
        prov_out_agri[colName] = prov_tmp.loc[:,'sum'] / pow(10,6)
    
    elif 'fores' in key:
        print(key)
        colName = key + '_Mt'
        prov_out_fores[colName] = prov_tmp.loc[:,'sum'] / pow(10,6)
        
    elif 'eCrops' in key:
        print(key)
        colName = key + '_Mt'
        prov_out_eCrops[colName] = prov_tmp.loc[:,'sum'] / pow(10,6)
ncfile.close()

# BioScenario
ncfile = nc.Dataset('../outputData/NetCDF/BioEnergy_scenario.nc','r',format='NETCDF4')
ncfile_keys = list(ncfile.variables.keys())
for key in ncfile_keys:
    tmp = np.array(ncfile.variables[key])
    prov_tmp = zonalByProvin_Table(chinaProv, chinaTable, tmp)
    colName = key + '_EJ' # MJ, GJ, TJ, EJ
    prov_out_bioScen[colName] = prov_tmp.loc[:,'sum'] / pow(10,9)
ncfile.close()

# write provincial results to excel file
writer = pd.ExcelWriter("../outputData/Excel/provin_stat_map.xlsx", 
                        engine="openpyxl", mode='w')
prov_out_land.to_excel(writer, sheet_name="margin_land")
prov_out_agri.to_excel(writer, sheet_name="agri_res")
prov_out_fores.to_excel(writer, sheet_name="fores_res")
prov_out_eCrops.to_excel(writer, sheet_name="eCrops")
prov_out_bioScen.to_excel(writer, sheet_name="bioScen")
writer.close()
del writer


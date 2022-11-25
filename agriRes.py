# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 09:48:11 2022

@author: vicke
"""

import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import netCDF4 as nc

# ----------- 1. calculate agricultural residues for each province ------------
sta = pd.read_excel('../inputData/agriSta.xlsx',sheet_name='1_Yield')
para_rpr = pd.read_excel('../inputData/agriSta.xlsx',sheet_name='2_para_RPR')
para_cr = pd.read_excel('../inputData/agriSta.xlsx',sheet_name='2_para_CR')
straw = sta.loc[:,'rice':'potatoes'] * para_rpr.loc[:,'rice':'potatoes'] * para_cr.loc[:,'rice':'potatoes']
straw = straw.set_index(sta['OBJECTID'])

# ------------ 2. allocate them on different land use type --------------
# 2.1 read map: China province map; land use map; npp map;
ds = gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
chinaProv = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/LandUseType/landUse2018.tif')
landUse = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/Npp/npp_2015.tif')
npp = ds.GetRasterBand(1).ReadAsArray()

# 2.2 paddy land (rice)
land_paddy = landUse.copy()
land_paddy[land_paddy!=11] = 0
land_paddy[land_paddy!=0] = 1

# (1) calculate for each province
straw_grid_paddy = np.zeros((4000,7000),dtype=float)
for i in range(1,35):
    # these provinces do not have straw data
    if i in [2, 28,30]:
        continue
    # extract each province
    prov_mask = chinaProv.copy()
    prov_mask[prov_mask!=i] = 0
    prov_mask[prov_mask!=0] = 1
    # statistic data
    prov_straw = straw.loc[i,'rice']
    # npp in {province, paddyland}
    npp_prov_paddy = prov_mask * land_paddy * npp
    npp_provSum = np.sum(npp_prov_paddy)
    # if this province does not have npp value, average the results
    if npp_provSum == 0: 
        npp_prov = prov_mask * npp
        straw_xy = prov_straw * npp_prov / np.sum(npp_prov)
        straw_grid_paddy = straw_grid_paddy + straw_xy
        continue
    straw_xy = prov_straw * npp_prov_paddy / npp_provSum
    straw_grid_paddy = straw_grid_paddy + straw_xy

# 2.3 dry land (other grains)
land_dry = landUse.copy()
land_dry[land_dry!=12] = 0
land_dry[land_dry!=0] = 1
straw_grid_dry = np.zeros((4000,7000),dtype=float)

# (1) calculate for each province
for i in range(1,35):
    if i in [2,28,30]:
        continue
    # extract each province
    prov_mask = chinaProv.copy()
    prov_mask[prov_mask!=i] = 0
    prov_mask[prov_mask!=0] = 1
    # statistic data
    prov_straw = np.sum(straw.loc[i,'maize':'potatoes'])
    # npp in {province, dryland}
    npp_prov_dry = prov_mask * land_dry * npp
    npp_provSum = np.sum(npp_prov_dry)
    if npp_provSum == 0:
        npp_prov = prov_mask * npp
        straw_xy = prov_straw * npp_prov / np.sum(npp_prov)
        straw_grid_dry = straw_grid_dry + straw_xy
        continue
    straw_xy = prov_straw * npp_prov_dry / npp_provSum
    straw_grid_dry = straw_grid_dry + straw_xy
    
# 2.4 merge two dataset  
straw_grid = straw_grid_dry + straw_grid_paddy

# 2.5 write results in nc file
ncfile = nc.Dataset('../outputData/bioenergy.nc','a',format='NETCDF4')

# create variable
nc_agri = ncfile.createVariable('agri_Res', np.float32, ('lat','lon'))
nc_agri.units = 'tonne/km2'
nc_agri.standard_name = 'agricultural_residues_production'
nc_agri[:,:] = straw_grid

# close file
ncfile.close()


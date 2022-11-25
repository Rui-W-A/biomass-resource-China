# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 15:08:11 2022

@author: vicke
"""

import numpy as np
import pandas as pd
from osgeo import gdal
import matplotlib.pyplot as plt
import netCDF4 as nc

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

# ------------ 1. calculate forestry residues for each province ----------------
forest_1 = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='1_forest',index_col='OBJECTID')
forest_2 = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='2_bamboo',index_col='OBJECTID')
forest_3 = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='3_intermediate_cutting',index_col='OBJECTID')
forest_4 = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='4_waste_wood_recycling',index_col='OBJECTID')
para = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='3_para',index_col='OBJECTID')
prov_ID =  forest_1.index.values

# 1.1 wood residues
merWood = forest_1.loc[:,'merWood_m3'] * 0.618 * (1-0.7717)/0.7717
fireWood = forest_1.loc[:,'fireWood_m3'] * 0.618
logWood = forest_1.loc[:,'log_m3'] * 0.618 * 0.6
res_wood = np.sum([merWood,fireWood,logWood], axis=0)

# 1.2 bamboo residues
bamboo = forest_2.loc[:,'bamboo'] * 0.015 * 0.3807 + forest_2.loc[:,'bamboo'] * 0.015 * 0.62
res_bamboo = np.array(bamboo)

# 1.3 intermediate cutting
fuelwood = forest_3.loc[:,'fuelForest'] * para.loc[:,'fuel1'] * para.loc[:,'fuel2']
timberwood = forest_3.loc[:,'timberForest'] * para.loc[:,'timber1'] * para.loc[:,'timber2']
protectwood = forest_3.loc[:,'protectForest']  * para.loc[:,'protect1'] * para.loc[:,'protect2']
economicwood = forest_3.loc[:,'economicForest']  * para.loc[:,'economic1'] * para.loc[:,'economic2']
specialwood = forest_3.loc[:,'specialForest'] * para.loc[:,'special1'] * para.loc[:,'special2']
streetwood = forest_3.loc[:,'streetTree'] * para.loc[:,'street1'] * para.loc[:,'street2']
res_intermediate = np.sum([fuelwood,timberwood,protectwood,economicwood,specialwood,streetwood], axis=0)

# 1.4 waste wood recycling
waste = forest_4.loc[:,'log_m3'] * 0.618 * 0.65
res_waste = np.array(waste)

# 1.5 merge all results
forestResi = pd.DataFrame([res_wood, res_bamboo, res_intermediate, res_waste],
                          columns = prov_ID,index=['wood','bamboo','intermediate','waste']).T

# --------------  2. spatial calculation and allocation --------------------
# ---- 2.1 read map: China province map; land use map; npp map;
ds = gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
chinaProv = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/LandUseType/landUse2018.tif')
landUse = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/Npp/npp_2015.tif')
npp = ds.GetRasterBand(1).ReadAsArray()

# ---- 2.2 intermediate cutting: res_intermediate on land type in {21,24}
land_21_24 = landUse.copy()
land_21_24[(land_21_24 == 21)|(land_21_24 == 24)] = 1
land_21_24[land_21_24!=1] = 0

# (1) calculate for each province
forestry_grid_inter1 = np.zeros((4000,7000),dtype=float)
for i in range(1,35):
    # these provinces do not have straw data
    if i in [2, 28,30]:
        continue
    # extract each province
    prov_mask = chinaProv.copy()
    prov_mask[prov_mask!=i] = 0
    prov_mask[prov_mask!=0] = 1
    # statistic data
    prov_fores = forestResi.loc[i,'intermediate']
    # npp in {province, paddyland}
    npp_prov_21_24 = prov_mask * land_21_24 * npp
    npp_provSum = np.sum(npp_prov_21_24)
    # if this province does not have npp value, average the results
    if npp_provSum == 0: 
        npp_prov = prov_mask * npp
        fores_xy = prov_fores * npp_prov / np.sum(npp_prov)
        forestry_grid_inter1 = forestry_grid_inter1 + fores_xy
        continue
    fores_xy = prov_fores * npp_prov_21_24 / npp_provSum
    forestry_grid_inter1 = forestry_grid_inter1 + fores_xy

# ---- 2.3 intermediate cutting: shrub on land = 22
land_22 = landUse.copy()
land_22[land_22 == 22] = 1
land_22[land_22!=1] = 0

# (1) calculate for each province
forestry_grid_inter2 = np.zeros((4000,7000),dtype=float)
for i in range(1,35):
    # these provinces do not have straw data
    if i in [2, 28,30]:
        continue
    # extract each province
    prov_mask = chinaProv.copy()
    prov_mask[prov_mask!=i] = 0
    prov_mask[prov_mask!=0] = 1
    # parameters
    para_1 = para.loc[i,'shrub1']
    para_2 = para.loc[i,'shrub2']
    
    # intermediate cutting for shrub in province i
    shrub_cutting_prov = prov_mask * land_22 * para_1 * para_2 * chinaArea
    forestry_grid_inter2 = forestry_grid_inter2 + shrub_cutting_prov

# ---- 2.4 intermediate cutting: sparse on land = 23
land_23 = landUse.copy()
land_23[land_23 == 23] = 1
land_23[land_23!=1] = 0

# (1) calculate for each province
forestry_grid_inter3 = np.zeros((4000,7000),dtype=float)
for i in range(1,35):
    # these provinces do not have straw data
    if i in [2, 28, 30]:
        continue
    # extract each province
    prov_mask = chinaProv.copy()
    prov_mask[prov_mask!=i] = 0
    prov_mask[prov_mask!=0] = 1
    # parameters
    para_1 = para.loc[i,'sparse1']
    para_2 = para.loc[i,'sparse2']
    
    # intermediate cutting for shrub in province i
    shrub_cutting_prov = prov_mask * land_23 * para_1 * para_2 * chinaArea
    forestry_grid_inter3 = forestry_grid_inter3 + shrub_cutting_prov
    
# ---- 2.5 others forestry residues on {21,22,23,24}
land_fores = landUse.copy()
land_fores[(land_fores == 21)|(land_fores == 22)|(land_fores == 23)|(land_fores == 24)] = 1
land_fores[land_fores!=1] = 0

# (1) calculate for each province
forestry_grid_others = np.zeros((4000,7000),dtype=float)
for i in range(1,35):
    # these provinces do not have straw data
    if i in [2,28,30]: # 2: Aomen
        continue
    # extract each province
    prov_mask = chinaProv.copy()
    prov_mask[prov_mask!=i] = 0
    prov_mask[prov_mask!=0] = 1
    # statistic data
    prov_fores = np.sum(forestResi.loc[i,['wood','bamboo','waste']])
    # npp in {province, paddyland}
    npp_prov_fores = prov_mask * land_fores * npp
    npp_provSum = np.sum(npp_prov_fores)
    # if this province does not have npp value, average the results
    if npp_provSum == 0: 
        npp_prov = prov_mask * npp
        fores_xy = prov_fores * npp_prov / np.sum(npp_prov)
        forestry_grid_others = forestry_grid_others + fores_xy
        continue
    fores_xy = prov_fores * npp_prov_fores / npp_provSum
    forestry_grid_others = forestry_grid_others + fores_xy

# ---------------  3. forestry residues -----------------
forestry_grid = forestry_grid_inter1 + forestry_grid_inter2 + forestry_grid_inter3 + forestry_grid_others

# write results in nc file
ncfile = nc.Dataset('../outputData/bioenergy.nc','a',format='NETCDF4')

# create variable
nc_agri = ncfile.createVariable('fores_Res', np.float32, ('lat','lon'))
nc_agri.units = 'tonne/km2'
nc_agri.standard_name = 'forestry_residues_production'
nc_agri[:,:] = forestry_grid

# close file
ncfile.close()



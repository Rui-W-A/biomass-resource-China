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
forest_3 = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='3_logging',index_col='OBJECTID')
# forest_4 = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='4_waste_wood_recycling',index_col='OBJECTID')
para = pd.read_excel('../inputData/foresSta.xlsx',sheet_name='3_para',index_col='OBJECTID')
prov_ID =  forest_1.index.values
ncfile = nc.Dataset('../outputData/NetCDF/BioRes_breakdown.nc','a',format='NETCDF4')
forestry_all = np.zeros((4000,7000),dtype=float)

# 1.1 wood residues
merWood = forest_1.loc[:,'merWood_m3'] * 0.618 * (1-0.7717)/0.7717
fireWood = forest_1.loc[:,'fireWood_m3'] * 0.618
logWood = forest_1.loc[:,'log_m3'] * 0.618 * 0.6
# res_wood = np.sum([merWood,fireWood,logWood], axis=0) 

# 1.2 bamboo residues
bamboo = forest_2.loc[:,'bamboo'] * 0.015 * 0.3807 + forest_2.loc[:,'bamboo'] * 0.015 * 0.62
res_bamboo = np.array(bamboo)

# 1.3 intermediate cutting
fuelwood = forest_3.loc[:,'fuelForest'] * para.loc[:,'fuel1'] * para.loc[:,'fuel2']
timberwood = forest_3.loc[:,'timberForest'] * para.loc[:,'timber1'] * para.loc[:,'timber2']
protectwood = forest_3.loc[:,'protectForest']  * para.loc[:,'protect1'] * para.loc[:,'protect2']
economicwood = forest_3.loc[:,'economicForest']  * para.loc[:,'economic1'] * para.loc[:,'economic2']
specialwood = forest_3.loc[:,'specialForest'] * para.loc[:,'special1'] * para.loc[:,'special2']
# streetwood = forest_3.loc[:,'streetTree'] * para.loc[:,'street1'] * para.loc[:,'street2']
# res_intermediate = np.sum([fuelwood,timberwood,protectwood,economicwood,specialwood,streetwood], axis=0)

# 1.4 waste wood recycling
# waste = forest_4.loc[:,'log_m3'] * 0.618 * 0.65
# res_waste = np.array(waste)

# 1.5 merge all results
# forestResi = pd.DataFrame([res_wood, res_bamboo, res_intermediate, res_waste],
#                           columns = prov_ID,index=['wood','bamboo','intermediate','waste']).T

# indexName = ['woodprocess_merchant', 'woodprocess_fire', 'woodprocess_log', 'woodprocess_bamboo',
#              'logging_fuel', 'logging_timber', 'logging_protect', 'logging_economic',
#              'logging_specialpurpose', 'logging_street']
# forestResi = pd.DataFrame([merWood, fireWood, logWood, bamboo, fuelwood,
#                            timberwood, protectwood, economicwood, specialwood, streetwood],
#                           columns = prov_ID, index = indexName).T

indexName = ['woodprocess_merchant', 'woodprocess_fire', 'woodprocess_log', 'woodprocess_bamboo',
             'logging_fuel', 'logging_timber', 'logging_protect', 'logging_economic',
             'logging_specialpurpose']
forestResi = pd.DataFrame([merWood, fireWood, logWood, bamboo, fuelwood,
                           timberwood, protectwood, economicwood, specialwood],
                          columns = prov_ID, index = indexName).T

# forestResi['fores_all'] = np.sum(forestResi, axis=1)

# save it in the excel file
chinaTable = pd.read_excel('../inputData/map/ChinaMask/ChinaTable.xlsx', index_col='OBJECTID')
foresRes_merge = pd.merge(chinaTable, forestResi, left_index=True, right_index=True, how='left')
foresRes_merge = foresRes_merge.fillna(0)

writer = pd.ExcelWriter("../outputData/Excel/provin_stat_account.xlsx", engine="openpyxl", mode='a')
foresRes_merge.to_excel(writer, sheet_name="fores_res")
writer.close()
del writer

# --------------  2. spatial calculation and allocation --------------------
# ---- 2.1 read map: China province map; land use map; npp map;
ds = gdal.Open('../inputData/map/ChinaMask/ChinaProvince.tif')
chinaProv = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/LandUseType/landUse2018.tif')
landUse = ds.GetRasterBand(1).ReadAsArray()
ds = gdal.Open('../inputData/map/Npp/npp_2015.tif')
npp = ds.GetRasterBand(1).ReadAsArray()

npp[npp>800] = 0
landUse_forest = landUse.copy()
landUse_forest[(landUse_forest == 21)|(landUse_forest == 22)|(landUse_forest == 23)|(landUse_forest == 24)] = 1
landUse_forest[(landUse_forest == 21)] = 1
landUse_forest[landUse_forest!=1] = 0

# ---- 2.2 intermediate cutting: res_intermediate on land type in {21}
# land_21_24 = landUse.copy()
# land_21_24[(land_21_24 == 21)|(land_21_24 == 24)] = 1
# land_21_24[land_21_24!=1] = 0
land_21 = landUse.copy()
land_21[land_21==21] = 1
land_21[land_21!=1] = 0

# (1) calculate for each province
# thinningName = indexName [4:]
for j in indexName:
    forestry_grid_inter1 = np.zeros((4000,7000),dtype=float)
    for i in range(1,35):
        # these provinces do not have straw data
        if i in [2, 28, 30]:
            continue
        # extract each province
        prov_mask = chinaProv.copy()
        prov_mask[prov_mask!=i] = 0
        prov_mask[prov_mask!=0] = 1
        # statistic data
        prov_fores = forestResi.loc[i,j]
        # npp in {province, paddyland}
        npp_prov_21 = prov_mask * land_21 * npp
        npp_provSum = np.sum(npp_prov_21)
        # if this province does not have npp value, average the results
        if npp_provSum == 0: 
            print(i)
            npp_prov = prov_mask * npp * landUse_forest
            fores_xy = prov_fores * npp_prov / np.sum(npp_prov)
            forestry_grid_inter1 = forestry_grid_inter1 + fores_xy
            continue
        fores_xy = prov_fores * npp_prov_21 / npp_provSum
        forestry_grid_inter1 = forestry_grid_inter1 + fores_xy
        
    name_j = 'fores_' + j
    print(name_j)
    nc_fores = ncfile.createVariable(name_j, np.float32, ('lat','lon'))
    nc_fores.units = 'tonne'
    nc_fores.standard_name = name_j
    nc_fores[:,:] = forestry_grid_inter1
    
    forestry_all = forestry_all + forestry_grid_inter1

# ---- 2.3 intermediate cutting: shrub on land = 22
land_22 = landUse.copy()
land_22[land_22 == 22] = 1
land_22[land_22!=1] = 0

# (1) calculate for each province
forestry_grid_inter2 = np.zeros((4000,7000),dtype=float)
for i in range(1,35):
    # these provinces do not have straw data
    if i in [2, 28, 30]:
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
    
nc_fores = ncfile.createVariable('fores_logging_shrub', np.float32, ('lat','lon'))
nc_fores.units = 'tonne'
nc_fores.standard_name = 'fores_logging_shrub'
nc_fores[:,:] = forestry_grid_inter2
forestry_all = forestry_all + forestry_grid_inter2

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
    
nc_fores = ncfile.createVariable('fores_logging_sparse', np.float32, ('lat','lon'))
nc_fores.units = 'tonne'
nc_fores.standard_name = 'fores_logging_sparse'
nc_fores[:,:] = forestry_grid_inter3
forestry_all = forestry_all + forestry_grid_inter3

# ---- 2.5 others forestry residues on {21,22,23,24}
# land_fores = landUse.copy()
# land_fores[(land_fores == 21)|(land_fores == 22)|(land_fores == 23)|(land_fores == 24)] = 1
# land_fores[land_fores!=1] = 0

# # (1) calculate for each province
# woodProcessName = indexName [:4]
# for j in woodProcessName: 
#     forestry_grid_others = np.zeros((4000,7000),dtype=float)
#     for i in range(1,35):
#         # these provinces do not have straw data
#         if i in [2,28,30]: # 2: Aomen
#             continue
#         # extract each province
#         prov_mask = chinaProv.copy()
#         prov_mask[prov_mask!=i] = 0
#         prov_mask[prov_mask!=0] = 1
#         # statistic data
#         prov_fores = forestResi.loc[i,j]
#         # npp in {province, paddyland}
#         npp_prov_fores = prov_mask * land_fores * npp
#         npp_provSum = np.sum(npp_prov_fores)
#         # if this province does not have npp value, average the results
#         if npp_provSum == 0: 
#             print(i)
#             npp_prov = prov_mask * npp * landUse_forest
#             fores_xy = prov_fores * npp_prov / np.sum(npp_prov)
#             forestry_grid_others = forestry_grid_others + fores_xy
#             continue
#         fores_xy = prov_fores * npp_prov_fores / npp_provSum
#         forestry_grid_others = forestry_grid_others + fores_xy
        
#     name_j = 'fores_' + j
#     print(name_j)
#     nc_fores = ncfile.createVariable(name_j, np.float32, ('lat','lon'))
#     nc_fores.units = 'tonne'
#     nc_fores.standard_name = name_j
#     nc_fores[:,:] = forestry_grid_others
    
    # forestry_all = forestry_all + forestry_grid_others

# ---------------  3. forestry residues -----------------
nc_fores = ncfile.createVariable('fores_all', np.float32, ('lat','lon'))
nc_fores.units = 'tonne'
nc_fores.standard_name = 'fores_all'
nc_fores[:,:] = forestry_all
ncfile.close()
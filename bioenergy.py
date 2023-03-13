# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:03:54 2022

@author: vicke
"""

import numpy as np
import pandas as pd
import netCDF4 as nc

# four scenar, typeSource（EJ）
ncfile_ori = nc.Dataset('../outputData/NetCDF/BioRes_breakdown.nc', 'r', format='NETCDF4')

agri = np.array(ncfile_ori.variables['agri_all']) * 14.7 # GJ/tonne
fores = np.array(ncfile_ori.variables['fores_all']) * 17.3 # GJ/tonne
eCrop_margin = np.array(ncfile_ori.variables['eCrops_best_margin']) * 16.3 # GJ/tonne
eCrop_abandon = np.array(ncfile_ori.variables['eCrops_best_abandon']) * 16.3 # GJ/tonne

# maximum results
gridNum = 4000 * 7000
sce_bioenergy = np.zeros([gridNum, 4])
sce_bioenergy[:, 0] = agri.reshape(gridNum)
sce_bioenergy[:, 1] = fores.reshape(gridNum)
sce_bioenergy[:, 2] = eCrop_abandon.reshape(gridNum) # abandon cropland
sce_bioenergy[:, 3] = eCrop_margin.reshape(gridNum)  # marginal land

# four scenarios
AF_min = np.nanmax(sce_bioenergy[:,:2], axis=1) * 0.11 # 11% of Agri & Forestry
AF_max = np.nanmax(sce_bioenergy[:,:2], axis=1)  # Agri & Forestry
AFE_abandon = np.nanmax(sce_bioenergy[:,:3], axis=1) # A+F+E on abandon cropland
AFE_margin = np.nanmax(sce_bioenergy[:,[0,1,3]], axis=1) # A+F+E on marginal land
AFE_max = np.nanmax(sce_bioenergy[:,:4], axis=1) # A+F+E on full exploitation

AF_min = AF_min.reshape(4000,7000)
AF_max = AF_max.reshape(4000,7000)
AFE_abandon = AFE_abandon.reshape(4000,7000)
AFE_margin = AFE_margin.reshape(4000,7000)
AFE_max = AFE_max.reshape(4000,7000)

# write into nc file
ncfile = nc.Dataset('../outputData/NetCDF/BioEnergy_scenario.nc','w',format='NETCDF4')
ncfile.createDimension('lat', 4000) # latitude axis
ncfile.createDimension('lon', 7000)  # longitude axis

nc_AF_min = ncfile.createVariable('s5_AF_min', np.float32, ('lat','lon'))
nc_AF_min.units = 'GJ/km2'
nc_AF_min.standard_name = '11% of agricultural and forestry residues'
nc_AF_min[:,:] = AF_min

nc_AF_max = ncfile.createVariable('s4_AF_max', np.float32, ('lat','lon'))
nc_AF_max.units = 'GJ/km2'
nc_AF_max.standard_name = '100% of agricultural and forestry residues'
nc_AF_max[:,:] = AF_max

nc_AFE_min = ncfile.createVariable('s3_AFE_abandon', np.float32, ('lat','lon'))
nc_AFE_min.units = 'GJ/km2'
nc_AFE_min.standard_name = 'agricultural, forestry residues, and energy crop on abandoned cropland'
nc_AFE_min[:,:] = AFE_abandon

nc_AFE_max = ncfile.createVariable('s2_AFE_margin', np.float32, ('lat','lon'))
nc_AFE_max.units = 'GJ/km2'
nc_AFE_max.standard_name = 'agricultural, forestry residues, and energy crop on marginal land'
nc_AFE_max[:,:] = AFE_margin

nc_AFE_max = ncfile.createVariable('s1_AFE_max', np.float32, ('lat','lon'))
nc_AFE_max.units = 'GJ/km2'
nc_AFE_max.standard_name = 'agricultural, forestry residues, and energy crop on all available marginal land'
nc_AFE_max[:,:] = AFE_max

ncfile.close()

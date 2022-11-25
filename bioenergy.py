# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 14:03:54 2022

@author: vicke
"""

import numpy as np
import pandas as pd
import netCDF4 as nc

# 四组情景，typeSource（EJ）
ncfile = nc.Dataset('../outputData/bioenergy.nc','a',format='NETCDF4')

agri = np.array(ncfile.variables['agri_Res']) * 14.7 # GJ/tonne
fores = np.array(ncfile.variables['fores_Res']) * 17.3 # GJ/tonne
eCrop_max = np.array(ncfile.variables['eCrop_margin']) * 16.3 # GJ/tonne
eCrop_min = np.array(ncfile.variables['eCrop_abandon']) * 16.3 # GJ/tonne

# maximum results
gridNum = 4000 * 7000
sce_bioenergy = np.zeros([gridNum, 8])
sce_bioenergy[:, 0] = agri.reshape(gridNum)
sce_bioenergy[:, 1] = fores.reshape(gridNum)
sce_bioenergy[:, 2] = eCrop_max.reshape(gridNum)
sce_bioenergy[:, 3] = eCrop_min.reshape(gridNum)

# four scenarios
AF_min = np.nanmax(sce_bioenergy[:,:2], axis=1) * 0.11 # 11% of Agri & Forestry
AF_max = np.nanmax(sce_bioenergy[:,:2], axis=1)  # Agri & Forestry
AFE_min = np.nanmax(sce_bioenergy[:,:3], axis=1) # A+F+E on abandon cropland
AFE_max = np.nanmax(sce_bioenergy[:,:4], axis=1) # A+F+E on full exploitation

AF_min = AF_min.reshape(4000,7000)
AF_max = AF_max.reshape(4000,7000)
AFE_min = AFE_min.reshape(4000,7000)
AFE_max = AFE_max.reshape(4000,7000)

# write into nc file
ncfile = nc.Dataset('../outputData/bioenergy.nc','a',format='NETCDF4')
nc_AF_min = ncfile.createVariable('AF_min', np.float32, ('lat','lon'))
nc_AF_min.units = 'GJ/km2'
nc_AF_min.standard_name = '11% of agricultural and forestry residues'
nc_AF_min[:,:] = AF_min

nc_AF_max = ncfile.createVariable('AF_max', np.float32, ('lat','lon'))
nc_AF_max.units = 'GJ/km2'
nc_AF_max.standard_name = '100% of agricultural and forestry residues'
nc_AF_max[:,:] = AF_max

nc_AFE_min = ncfile.createVariable('AFE_min', np.float32, ('lat','lon'))
nc_AFE_min.units = 'GJ/km2'
nc_AFE_min.standard_name = 'agricultural, forestry residues, and energy crop on abandon cropland'
nc_AFE_min[:,:] = AFE_min

nc_AFE_max = ncfile.createVariable('AFE_max', np.float32, ('lat','lon'))
nc_AFE_max.units = 'GJ/km2'
nc_AFE_max.standard_name = 'agricultural, forestry residues, and maximum energy crop'
nc_AFE_max[:,:] = AFE_max

ncfile.close()

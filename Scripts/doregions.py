#!/usr/bin/env python
import numpy as np
import pandas as pd
import xarray as xr
import scipy.interpolate as si

regions = pd.read_csv('/Users/meilers/MITinternship/Data/Regions/regions.csv', header=None).values
reg_lon = np.loadtxt('/Users/meilers/MITinternship/Data/Regions/lons.csv', delimiter=',')
reg_lat = np.loadtxt('/Users/meilers/MITinternship/Data/Regions/lats.csv', delimiter=',')

# create inteprolator
I = si.RegularGridInterpolator((reg_lat, reg_lon), regions, 'nearest',
                               bounds_error=False, fill_value=None)

diazotroph_observations = pd.read_csv(r'/Users/meilers/MITinternship/Data/Tang_and_Cassar-2019/nifH_Gene_Integral_mod.csv')
#print(diazotroph_observations)

nifH_database = pd.DataFrame(diazotroph_observations, columns = ['LONGITUDE',
                                                                 'LATITUDE',
                                                                 'YEAR',
                                                                 'MONTH',
                                                                 'Trichodesmium nifH Gene (x106 copies m-2)',
                                                                 'UCYN-A nifH Gene (x106 copies m-2)',
                                                                 'UCYN-B nifH Gene (x106 copies m-2)',
                                                                 'UCYN-C nifH Gene (x106 copies m-2)',
                                                                 'Richelia nifH Gene (x106 copies m-2)',
                                                                 'Calothrix nifH Gene  (x106 copies m-2)',
                                                                 'Gamma nifH Gene (x106 copies/m3)'])

lon_nifH = diazotroph_observations['LONGITUDE'].astype(np.float32)
lat_nifH = diazotroph_observations['LATITUDE'].astype(np.float32)

n_d = len(lon_nifH)
latlon = np.zeros((n_d, 2))
latlon[:,0] = lat_nifH
# reg_lon has range [-180,180)
latlon[:,1] = np.mod(lon_nifH + 180, 360.) - 180

# interpolate
reg_nifH = I(latlon)


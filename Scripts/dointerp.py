#!/usr/bin/env python
import numpy as np
import xarray as xr
import scipy.interpolate as si

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')
lon = grid.X.values
lat = grid.Y.values
vals = grid.XC.values  # for testing - replace by data to be interpolated

# create inteprolator
I = si.RegularGridInterpolator((lat, lon), vals, 'nearest',
                               bounds_error=False, fill_value=None)

# where to evaluate data
lon_d = np.r_[0., 10., 20.]
lat_d = np.r_[-80., 0., 10.]

n_d = len(lon_d)
latlon = np.zeros((n_d, 2))
latlon[:,0] = lat_d
latlon[:,1] = np.mod(lon_d, 360.)

# interpolate
vals_d = I(latlon)



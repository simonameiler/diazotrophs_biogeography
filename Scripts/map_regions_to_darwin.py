import numpy as np
import xarray as xr
import pandas as pd
import scipy.interpolate as si

grid = xr.open_dataset('grid.nc')
lon = grid.X
lat = grid.Y

lm = grid.HFacC[0].values == 0

regions = pd.read_csv('regions.csv', header=None).values
#reg_lon = np.loadtxt('lons.csv', delimiter=',').astype(int)
#reg_lat = np.loadtxt('lats.csv', delimiter=',').astype(int)

reg_lat = np.arange(-89., 90., 2.)
reg_lon = np.arange(-181., 178., 2.)

# extend for periodicity in lon
reg_lon_extended = np.r_[reg_lon-360, reg_lon, reg_lon+360]
regions_extended = np.c_[regions, regions, regions]

# make 2d
yr, xr = np.meshgrid(reg_lat, reg_lon_extended, indexing='ij')

# do not use land values
w = regions_extended != 0

# make target coordinates 2d
xo, yo = np.meshgrid(lat, (lon+182)%360-182, indexing='ij')

# map
reg_darwin = si.griddata((yr[w],xr[w]), regions_extended[w], (xo.ravel(), yo.ravel()), 'nearest')
reg_darwin = reg_darwin.reshape(160, 360)
reg_darwin[lm] = 0

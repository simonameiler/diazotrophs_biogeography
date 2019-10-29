# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%% Import Python modules
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%% Load data
data = xr.open_dataset('/Users/meilers/MITinternship/Data/scenario3.nc')

lon = data.lon
lat = data.lat

iron = data.scenarioIV
iron_mean = np.mean(iron, axis=0)
iron_unit = np.log(iron_mean[:,:])#*31536000000)

#%% variables and constants
f_atm = data.dst_emis_solfe
k = 0.1
Fe0 = 1e-5
N0 = 1.6
rpFe = 6.25e-5

#%% Plot
fig,ax = plt.subplots(figsize=(9,6))
c = ax.contourf(lon,lat,iron_mean,levels=np.linspace(0,1.8e-12))
plt.show()

#%% Plot nice map
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)}, sharey=True, figsize=(9,4))
#ax.imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
ax.coastlines(color='#888888',linewidth=1.5)
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
c0 = ax.contourf(lon,lat,iron_unit)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.grid()
fig.colorbar(c0)
#fig.subplots_adjust(right=0.85)
#cb = fig.add_axes([0.9,0.25,0.01,0.5])
#fig.colorbar(c0,cax=cb,label='Depth (m)')
#cb.set_label('[m]',rotation=90, position=(0.5,1))
plt.show()

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from netCDF4 import Dataset

ds = Dataset('/Users/meilers/MITinternship/Data/MarEDat20130403Diazotrophs.nc', 'r')

lon = ds.variables['LONGITUDE']
lat = ds.variables['LATITUDE']

obs = ds.variables['OBSERVATIONS']
abund = ds.variables['ABUNDANCE']

obs_int = np.sum(obs,axis=1)
abund_int = np.sum(abund,axis=1)

obs_tot = np.sum(obs_int,axis=0)
abund_tot = np.sum(abund_int,axis=0)

#%% Create presence/absence matrix

pres_abs = np.zeros(180,360,1)

for i in len(obs_tot[0]):
    for j in len(obs_tot[:,0]):
        if obs_tot[i,j] != 0:
            pres_abs[i,j][2] == 1
        else:
            pres_abs[i,j][2] == 0

#%%
        
mask = np.zeros_like(obs_tot)
mask[obs_tot>0] = 1

find = np.where(mask==1)

#%%
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3)) # Map projection
ax.coastlines(color='k',linewidth=1)        # Adds coastline to map at highest resolution
ax.scatter(lon, lat, c=obs_tot, cmap='viridis', s=abund_tot, linewidth=0, alpha=0.5,transfor=ccrs.PlateCarree())
#plt.scatter(lngArr, latArr, s=area, c=satLngArr, alpha=0.5, transform=ccrs.Geodetic())                # Plot
plt.show()

#plt.scatter(lon, lat, label=None, c=np.sum(obs,axis=1), cmap='viridis', s=np.sum(abund,axis=1), linewidth=0, alpha=0.5)

#%%
obs.plot(kind="scatter",x="lon",y="lat",alpha=0.5)
plt.show()

#%% 
plt.scatter(obs_int[0],obs_int[1])

#%%
colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.plot(c)#,cmap=colmap)#,levels=levs[nut],cmap=colmap,extend='both')
#con = ax.contourf(non_z)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#cbar = plt.colorbar(c,ax=ax)
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()
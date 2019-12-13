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

#%% Load diazotroph data
ds = Dataset('/Users/meilers/MITinternship/Data/MarEDat20130403Diazotrophs.nc', 'r')

# extract variables which are needed and convert/integrate
lon = ds.variables['LONGITUDE']
lat = ds.variables['LATITUDE']

obs = ds.variables['OBSERVATIONS']
abund = ds.variables['ABUNDANCE']

obs_int = np.sum(obs[:,0:6,:,:],axis=1)
abund_int = np.sum(abund[:,0:6,:,:],axis=1)

obs_tot = np.sum(obs_int,axis=0)
abund_tot = np.sum(abund_int,axis=0)

#%% Load nutrient data
# iron dust flux
fein = np.fromfile('/Users/meilers/MITinternship/Data/mahowald2009_solubile_current_smooth_oce_mth-2d.bin', '>f4').reshape(12,160,360)

# model simulation output
Nutr = xr.open_dataset('/Users/meilers/MITinternship/Data/Nutr_tend.0000014400.nc')

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/supply50m.nc')

lon1 = grid.lon   #needed for plotting
lat1 = grid.lat    #needed for plotting

#%% 
# transport terms
NO3 = Nutr['gTr04']
PO4 = Nutr['gTr05']
FeT = Nutr['gTr07']

# other terms: remineralization and dust in case of Fe
NO3_o = Nutr['gGUD04']
PO4_o = Nutr['gGUD05']
FeT_o = Nutr['gGUD07']

sec = 31557600  #seconds per year (365.25*86400)
depth = 6 #RF = 0, -10, -20, -35, -55, -75, -100, -135, -185, (...)

NO3_trans = np.sum(NO3[0,0:depth,:,:],axis=0)*sec
PO4_trans = np.sum(PO4[0,0:depth,:,:],axis=0)*sec
FeT_trans = np.sum(FeT[0,0:depth,:,:],axis=0)*sec

NO3_other = np.sum(NO3_o[0,0:depth,:,:],axis=0)*sec
PO4_other = np.sum(PO4_o[0,0:depth,:,:],axis=0)*sec
FeT_other = np.sum(FeT_o[0,0:depth,:,:],axis=0)*sec

f_dust = np.sum(fein,axis=0)*sec

N = np.add(NO3_trans,NO3_other)
P = np.add(PO4_trans,PO4_other)
Fe = np.add(FeT_trans,FeT_other)#,f_dust) #including dust

#%% Calculate ratios

# Bioavailable nutrient supply --> bio_PN, bio_FeN
# Constants
rpP = 0.0625
rpFe = 6.25e-5
k = 0.1

#bio_PN = (np.divide(P,N))*(1/rpP)
#bio_FeN = (np.divide(Fe,N))*(1/rpFe)

bio_PN = (np.divide(PO4_trans,NO3_trans))*(1/rpP)
F = np.add((FeT_trans*k),f_dust)
bio_FeN = (np.divide(F,(NO3_trans*k)))*(1/rpFe)

#%% Create presence/absence matrix

#pres_abs = np.zeros(180,360,1)

#for i in len(obs_tot[0]):
#    for j in len(obs_tot[:,0]):
#        if obs_tot[i,j] != 0:
#            pres_abs[i,j][2] == 1
#        else:
#            pres_abs[i,j][2] == 0

#%% Create a mask: 0 = absence, 1 = presence
        
mask1 = np.zeros_like(obs_tot)
mask1[obs_tot>0] = 1

mask2 = np.zeros_like(abund_tot)
mask2[abund_tot>0] = 1

find1 = np.where(mask1==1)
find2 = np.where(mask2==1)

plt.plot(find1[1],find1[0],'.',color='b')
plt.plot(find2[1],find2[0],'.',color='r')
plt.show()

#%%
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3)) # Map projection
ax.coastlines(color='k',linewidth=1)        # Adds coastline to map at highest resolution
ax.scatter(lon, lat, c=obs_tot, cmap='viridis', s=abund_tot, linewidth=0, alpha=0.5,transfor=ccrs.PlateCarree())
#plt.scatter(lngArr, latArr, s=area, c=satLngArr, alpha=0.5, transform=ccrs.Geodetic())                # Plot
plt.show()

#plt.scatter(lon, lat, label=None, c=np.sum(obs,axis=1), cmap='viridis', s=np.sum(abund,axis=1), linewidth=0, alpha=0.5)

#%%
nut = 1    #chose nutrient here: 0=Fe, 1=N, 2=P
nutrient = [bio_PN,bio_FeN]
label_nut = ['N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)','Fe (mol/m$^{2}$/y)']
name_nut = ['P:N','Fe:N']

levs_PN = np.linspace(0,1000,11)
levs_FeN = np.linspace(-3,3,13)
levs = [levs_PN,levs_FeN]

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(15,5))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon1,lat1,np.log(nutrient[nut]),levels=levs[nut],cmap=colmap,extend='both')
#con1 = ax.contour(lon1,lat1,nutrient[nut],levels=[1.2],colors='k',linewidths=1,linstyle='solid')
#con2 = ax.contour(lon1,lat1,nutrient[nut],levels=[2.5],colors='r',linewidths=1,linstyle='solid')
plt.plot(lon[find1[1]],lat[find1[0]],'.',color='b')
plt.plot(lon[find2[1]],lat[find2[0]],'.',color='g')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label(''+str(name_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()
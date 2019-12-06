# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%% Load data

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/supply50m.nc')
area_info = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')

lon = grid.lon   #needed for plotting
lat = grid.lat    #needed for plotting

area = area_info.rA

# Read in diazotroph data from model - same setup as for nutrients
# 5 diazotroph species. according to Steph TRAC30 to TRAC34

months_vec = range(0,12)
#months_list = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mon_list = ['26160','26400','26640','26880','27120','27360','27600','27840','28080','28320','28560','28800']

diaz1 = np.zeros((len(months_vec),23,160,360))
diaz2 = np.zeros((len(months_vec),23,160,360))
diaz3 = np.zeros((len(months_vec),23,160,360))
diaz4 = np.zeros((len(months_vec),23,160,360))
diaz5 = np.zeros((len(months_vec),23,160,360))

# Open dataset
i = 0
for month in months_vec:
    diaz = xr.open_dataset('/Users/meilers/MITinternship/Data/run19_33_MONTH/3d.00000'+str(mon_list[month])+'.nc')
    diaz1[month,:,:,:] = diaz.TRAC30.values[:,:,:,:]
    diaz2[month,:,:,:] = diaz.TRAC31.values[:,:,:,:]
    diaz3[month,:,:,:] = diaz.TRAC32.values[:,:,:,:]
    diaz4[month,:,:,:] = diaz.TRAC33.values[:,:,:,:]
    diaz5[month,:,:,:] = diaz.TRAC34.values[:,:,:,:]
    diaz.close()
    print('Month: '+str(month))
    i += 1
# Sum up diazotroph data into one array
diaz = diaz1 + diaz2 + diaz3 + diaz4 + diaz5

#%% MAREDAT - Load and condense diazotroph data
ds = xr.open_dataset('/Users/meilers/MITinternship/Data/MarEDat20130403Diazotrophs.nc',decode_times=False)

# extract variables which are needed and convert/integrate
lon_d = ds['LONGITUDE']
lat_d = ds['LATITUDE'][10:-10] # to match the latitude of the nutrient data

obs = ds['OBSERVATIONS']
abund = ds['ABUNDANCE']
bm = ds['BIOMASS']
nifH = ds['nifHbiom']

obs_tot = np.sum(obs[:,0:6,10:-10,:],axis=(0,1))
abund_tot = np.sum(abund[:,0:6,10:-10,:],axis=(0,1))
bm_tot = np.sum(bm[:,:,:,:],axis=(0,1))
nifH_tot = np.sum(nifH[:,0:6,10:-10,:],axis=(0,1))


#%% define constants and dz
sec = 31557600  #seconds per year (365.25*86400)
dz = [10,10,15,20,20,25] #,35,50,75,100] #(...) dz between two depth layers
depth = [0,10,20,35,55,75,100,135,185,260,360,510,710,985,1335,1750,2200,2700,3200,3700,4200,4700,5200,5700]

# little loop to calculate all the dz's
dz_all = np.zeros(len(depth)-1)
for i in range(len(dz_all)):
    dz_all[i] = depth[i+1]-depth[i]
    
#%% sum nutrients up over depth and multiply with corresponding dz and sec

diaz_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz_int[:,i,:,:] = diaz[:,i,:,:]*dz_all[i]*sec  
    print(np.max(diaz_int[:,i,:,:]))
    print(i)
diaz_int = np.sum(diaz_int,axis=(0,1))

diaz_int_100 = np.zeros((12,6,160,360))
for i in range(len(dz)):
    diaz_int_100[:,i,:,:] = diaz[:,i,:,:]*dz[i]*sec  
    print(np.max(diaz_int_100[:,i,:,:]))
    print(i)
diaz_int_100 = np.sum(diaz_int_100,axis=(0,1))

#%% mask where P:N OR Fe:N is sufficient to support diazotrophs
mask = np.where((diaz_int > 31557600), 1, 0)
mask_out = np.where((diaz_int < 31557600), 1, 0)

#%% Manipulate diazotroph data
     
diaz_obs = np.zeros_like(obs_tot)
diaz_obs[obs_tot>0] = 1

diaz_abund = np.zeros_like(abund_tot)
diaz_abund[abund_tot>0] = 1

diaz_bm = np.zeros_like(bm_tot)
diaz_bm[bm_tot>0] = 1

diaz_nifH = np.zeros_like(nifH_tot)
diaz_nifH[nifH_tot>0] = 1

find_obs = np.where(diaz_obs==1)
find_abund = np.where(diaz_abund==1)
find_bm = np.where(diaz_bm==1)
find_nifH = np.where(diaz_nifH==1)

absent_obs = np.where(diaz_obs-diaz_abund==1)

#pack the masks for the different diazotroph variables into one list
diaz_data_list = [find_obs,find_abund,find_bm,find_nifH,absent_obs]

#%% Quantify how many of the diazotrophs abundances are in the predicted province
# careful: make sure to get lon/lat of nutrients and diazotrophs consistent!!!
# correct the two scales of latitude to match one another. (lon would be the same but to avoid confusion
# I converted it too.) All we care about here is getting the right indices matching the lon, lat of both,
# diazotroph and nutrient data. 
list_idx = 1 #to chose which data from diaz_data_list to plot

lat_corr = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
lon_corr = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
IN = np.sum(mask[lat_corr,lon_corr])/len(lat_corr)
print(IN)

#%% Calculate accuracz for absences
list_idx = 4
lat_corr_abs = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
lon_corr_abs = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
OUT = np.sum(mask_out[lat_corr_abs,lon_corr_abs])/len(lat_corr_abs)
print(OUT)
  
#%% Plot diazotroph biomass simulated in Darwin (integrated over top 100m, and over entire depth range; summed up over 1 year)

col = plt.get_cmap('RdBu_r')

depth_lab = ['top 100m','entire depth range']
fig,ax = plt.subplots(2,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,8),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,2):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.15,0.9,''+(str(depth_lab[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax[0].contourf(lon,lat,diaz_int_100,levels=np.linspace(0,1.5e10,20),cmap=col)#,extend='both')
c1 = ax[1].contourf(lon,lat,diaz_int,levels=np.linspace(0,1.5e10,20),cmap=col)#,extend='both')
#plt.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='g')
#plt.plot(lon_d[find_bm[1]],lat_d[find_bm[0]],'.',color='r')
#plt.plot(lon_d[find_nifH[1]],lat_d[find_nifH[0]],'.',color='c')
#plt.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='m')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nodata_alldepth.png', bbox_inches='tight', dpi=300)

#%% Plot diazotroph biomass simulated in Darwin (integrated over top 100m, summed up over 1 year)

col = plt.get_cmap('RdBu_r')

depth_lab = ['top 100m','entire depth range']
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
ax.text(0.15,0.9,''+(str(depth_lab[0])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,diaz_int_100,levels=np.linspace(0,1.5e10,20),cmap=col)#,extend='both')
#plt.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='g')
#plt.plot(lon_d[find_bm[1]],lat_d[find_bm[0]],'.',color='r')
#plt.plot(lon_d[find_nifH[1]],lat_d[find_nifH[0]],'.',color='c')
#plt.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='m')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nodata_100m.png', bbox_inches='tight', dpi=300)

#%% Plot diazotroph biomass simulated in Darwin (integrated over top 100m, summed up over 1 year)

col = plt.get_cmap('RdBu_r')

depth_lab = ['top 100m','entire depth range']
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
ax.text(0.15,0.9,''+(str(depth_lab[0])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,diaz_int_100,levels=np.linspace(0,1.5e10,20),cmap=col)#,extend='both')
ax.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='orange',label='presence')
ax.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='m',label='absence')
ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_withdata_100m.png', bbox_inches='tight', dpi=300)


#%% Just for some quick plots
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(12,4))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon_d,lat_d,bm_tot,levels=np.linspace(0,1,21),cmap=cm.cm.haline,extend='both')
#plt.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='b')
#plt.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='g')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
#cbar.set_label(''+str(name_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()
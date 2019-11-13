# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%% Load data

# iron dust flux
fein = np.fromfile('/Users/meilers/MITinternship/Data/mahowald2009_solubile_current_smooth_oce_mth-2d.bin', '>f4').reshape(12,160,360)

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/supply50m.nc')
area_info = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')

lon = grid.lon   #needed for plotting
lat = grid.lat    #needed for plotting

area = area_info.rA

#%% Read in monthly nutrient data

months_vec = range(0,12)
months_list = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mon_list = ['26160','26400','26640','26880','27120','27360','27600','27840','28080','28320','28560','28800']

#timevector = np.zeros((1,1))
NH4_timeseries = np.zeros((len(months_vec),6,160,360))
NO2_timeseries = np.zeros((len(months_vec),6,160,360))
NO3_timeseries = np.zeros((len(months_vec),6,160,360))
PO4_timeseries = np.zeros((len(months_vec),6,160,360))
FeT_timeseries = np.zeros((len(months_vec),6,160,360))
S_DIN_timeseries = np.zeros((len(months_vec),6,160,360))
S_PO4_timeseries = np.zeros((len(months_vec),6,160,360))
S_Fe_timeseries = np.zeros((len(months_vec),6,160,360))

# Open dataset
i = 0
for month in months_vec:
    data = xr.open_dataset('/Users/meilers/MITinternship/Data/Nutrients_monthly/nutr_tend.'+str(months_list[month])+'.nc')
    NH4_timeseries[month,:,:,:] = data.gDAR04.values[:,0:6,:,:]
    NO2_timeseries[month,:,:,:] = data.gDAR03.values[:,0:6,:,:]
    NO3_timeseries[month,:,:,:] = data.gDAR02.values[:,0:6,:,:]
    PO4_timeseries[month,:,:,:] = data.gDAR05.values[:,0:6,:,:]
    FeT_timeseries[month,:,:,:] = data.gDAR06.values[:,0:6,:,:]
    S_DIN_timeseries[month,:,:,:] = data.S_DIN.values[:,0:6,:,:]
    S_PO4_timeseries[month,:,:,:] = data.S_PO4.values[:,0:6,:,:]
    S_Fe_timeseries[month,:,:,:] = data.S_Fe.values[:,0:6,:,:]
    data.close()
    print('Month: '+str(month))
    i += 1


#%% Read in diazotroph data from model

diaz1_timeseries = np.zeros((len(months_vec),6,160,360))
diaz2_timeseries = np.zeros((len(months_vec),6,160,360))
diaz3_timeseries = np.zeros((len(months_vec),6,160,360))
diaz4_timeseries = np.zeros((len(months_vec),6,160,360))
diaz5_timeseries = np.zeros((len(months_vec),6,160,360))

# Open dataset
i = 0
for month in months_vec:
    diaz = xr.open_dataset('/Users/meilers/MITinternship/Data/Diaz_mon/3d.00000'+str(mon_list[month])+'.nc')
    diaz1_timeseries[month,:,:,:] = diaz.TRAC30.values[:,0:6,:,:]
    diaz2_timeseries[month,:,:,:] = diaz.TRAC31.values[:,0:6,:,:]
    diaz3_timeseries[month,:,:,:] = diaz.TRAC32.values[:,0:6,:,:]
    diaz4_timeseries[month,:,:,:] = diaz.TRAC33.values[:,0:6,:,:]
    diaz5_timeseries[month,:,:,:] = diaz.TRAC34.values[:,0:6,:,:]
    diaz.close()
    print('Month: '+str(month))
    i += 1
    

#%% Load and condense nutrient data
# transport terms
NH4 = NH4_timeseries
NO2 = NO2_timeseries
NO3 = NO3_timeseries
PO4 = PO4_timeseries
FeT = FeT_timeseries

# other terms: remineralization and dust in case of Fe
S_DIN = S_DIN_timeseries
S_PO4 = S_PO4_timeseries
S_Fe = S_Fe_timeseries


# diazotroph data
diaz = diaz1_timeseries + diaz2_timeseries + diaz3_timeseries + diaz4_timeseries + diaz5_timeseries

# define constants and dz
sec = 1  #seconds per year (365.25*86400)
dz = [10,10,15,20,20,25] #,35,50,75,100] #(...) dz between two depth layers

#%% sum nutrients up over depth and multiply with corresponding dz and sec

NH4_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    NH4_int[:,i,:,:] = NH4[:,i,:,:]*dz[i]*sec
#    print(np.max(NH4_int[:,i,:,:]))
#    print(i)
NH4_int = np.sum(NH4_int,axis=1)

NO2_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    NO2_int[:,i,:,:] = NO2[:,i,:,:]*dz[i]*sec    
NO2_int = np.sum(NO2_int,axis=1)
 
NO3_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    NO3_int[:,i,:,:] = NO3[:,i,:,:]*dz[i]*sec    
NO3_int = np.sum(NO3_int,axis=1)

PO4_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    PO4_int[:,i,:,:] = PO4[:,i,:,:]*dz[i]*sec    
PO4_int = np.sum(PO4_int,axis=1)

FeT_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    FeT_int[:,i,:,:] = FeT[:,i,:,:]*dz[i]*sec    
FeT_int = np.sum(FeT_int,axis=1)

S_DIN_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    S_DIN_int[:,i,:,:] = S_DIN[:,i,:,:]*dz[i]*sec    
S_DIN_int = np.sum(S_DIN_int,axis=1)

S_PO4_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    S_PO4_int[:,i,:,:] = S_PO4[:,i,:,:]*dz[i]*sec    
S_PO4_int = np.sum(S_PO4_int,axis=1)

S_Fe_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    S_Fe_int[:,i,:,:] = S_Fe[:,i,:,:]*dz[i]*sec    
S_Fe_int = np.sum(S_Fe_int,axis=1)

diaz_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    diaz_int[:,i,:,:] = diaz[:,i,:,:]*dz[i]*sec    
diaz_int = np.sum(diaz_int,axis=1)

#%% Add up the different N species
N_int = NH4_int + NO3_int + NO2_int
N_int_o = S_DIN_int

#%% Set all negative values to zero

#from copy import deepcopy
#N_trans = deepcopy(N_int)
#P_trans = deepcopy(PO4_int)
#Fe_trans = deepcopy(FeT_int)
#N_remin = deepcopy(N_int_o)
#P_remin = deepcopy(S_PO4_int)
#Fe_other = deepcopy(S_Fe_int)
  
#%% Define transport, remin and total nutrient fluxes; add iron flux

N_tot = np.add(N_int,N_int_o)
P_tot = np.add(PO4_int,S_PO4_int)
Fe_tot = np.add(FeT_int,S_Fe_int)

#%% Calculate seasonal means of nutrients here

N_tot_MAM = np.mean(N_tot[2:4,:,:],axis=0)
P_tot_MAM = np.mean(P_tot[2:4,:,:],axis=0)
Fe_tot_MAM = np.mean(Fe_tot[2:4,:,:],axis=0)

#%% Calculate ratios

# Bioavailable nutrient supply --> bio_PN, bio_FeN
# Constants
rpP = 0.0625
rpFe = 6.25e-5
k = 0.1
ref_PN = 1.04
ref_FeN = 1.2

bio_PN_tot = (np.divide(P_tot[:,:,:],N_tot[:,:,:]))*(1/rpP)
bio_FeN_tot = (np.divide(Fe_tot[:,:,:],N_tot[:,:,:]))*(1/rpFe)

#%% then calculate the nutrient ratios with the seasonal means of March, April, May

bio_PN_MAM = (np.divide(P_tot_MAM[:,:],N_tot_MAM[:,:]))*(1/rpP)
bio_FeN_MAM = (np.divide(Fe_tot_MAM[:,:],N_tot_MAM[:,:]))*(1/rpFe)
#bio_PN = (P_remin/N_remin)*(1/rpP)
#bio_FeN = (Fe_other/N_remin)*(1/rpFe)

#%% Control plot nutrient fields - N_tot, P_tot, Fe_tot for 1 month

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(3,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,9),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,3):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cN = ax[0].contourf(lon,lat,N_tot[0,:,:],levels=np.linspace(-2e-04,2e-04,101),cmap=colmap,extend='both')
cP = ax[1].contourf(lon,lat,P_tot[0,:,:],levels=np.linspace(-5e-06,5e-06,101),cmap=colmap,extend='both')
cFe = ax[2].contourf(lon,lat,Fe_tot[0,:,:],levels=np.linspace(-6e-09,6e-09,101),cmap=colmap,extend='both')
#con = ax.contour(lon,lat,mask,color='r')
cbarN = plt.colorbar(cN,ax=ax[0])
cbarN = plt.colorbar(cP,ax=ax[1])
cbarN = plt.colorbar(cFe,ax=ax[2])
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Control plot nutrient ratios P:N

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,bio_PN_tot[0,:,:],levels=np.linspace(0.5,1.5,11),cmap=colmap,extend='both')
#con = ax.contour(lon,lat,mask,color='r')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Control plot nutrient ratios Fe:N

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,bio_FeN_tot[0,:,:],levels=np.linspace(0.5,1.5,11),cmap=colmap,extend='both')
#con = ax.contour(lon,lat,mask,color='r')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Control plot nutrient ratios P:N for March, April, May

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,bio_PN_MAM[:,:],levels=np.linspace(0.5,1.5,11),cmap=colmap,extend='both')
#con = ax.contour(lon,lat,mask,color='r')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Control plot nutrient ratios Fe:N for March, April, May

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,bio_FeN_MAM[:,:],levels=np.linspace(0.5,1.5,11),cmap=colmap,extend='both')
#con = ax.contour(lon,lat,mask,color='r')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()
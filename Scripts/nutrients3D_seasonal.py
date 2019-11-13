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
    data = xr.open_dataset('/Users/meilers/MITinternship/Data/Nutr_mon/nutr_tend.00000'+str(mon_list[month])+'.nc')
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
    
#%% Load and condense diazotroph data
ds = xr.open_dataset('/Users/meilers/MITinternship/Data/MarEDat20130403Diazotrophs.nc',decode_times=False)

# extract variables which are needed and convert/integrate
lon_d = ds['LONGITUDE']
lat_d = ds['LATITUDE'][10:-10] # to match the latitude of the nutrient data

obs = ds['OBSERVATIONS']
abund = ds['ABUNDANCE']
bm = ds['BIOMASS']
nifH = ds['nifHbiom']
nz_obs = ds['NON_ZERO_OBS']
nz_abund = ds['NON_ZERO_ABUND']
nz_bm = ds['NON_ZERO_BIOM']
nz_nifH = ds['NON_ZERO_nifH']

obs_tot = np.sum(obs[:,0:6,10:-10,:],axis=(0,1))
abund_tot = np.sum(abund[:,0:6,10:-10,:],axis=(0,1))
bm_tot = np.sum(bm[:,0:6,10:-10,:],axis=(0,1))
nifH_tot = np.sum(nifH[:,0:6,10:-10,:],axis=(0,1))
nz_obs_tot = np.sum(nz_obs[:,0:6,10:-10,:],axis=(0,1))
nz_abund_tot = np.sum(nz_abund[:,0:6,10:-10,:],axis=(0,1))
nz_bm_tot = np.sum(nz_bm[:,0:6,10:-10,:],axis=(0,1))
nz_nifH_tot = np.sum(nz_nifH[:,0:6,10:-10,:],axis=(0,1))

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

from copy import deepcopy
N_trans = deepcopy(N_int)
P_trans = deepcopy(PO4_int)
Fe_trans = deepcopy(FeT_int)
N_remin = deepcopy(N_int_o)
P_remin = deepcopy(S_PO4_int)
Fe_other = deepcopy(S_Fe_int)
  
#%% Define transport, remin and total nutrient fluxes; add iron flux

N_tot = np.add(N_trans,N_remin)
P_tot = np.add(P_trans,P_remin)
Fe_tot = np.add(Fe_trans,Fe_other)

#%% Calculate ratios

# Bioavailable nutrient supply --> bio_PN, bio_FeN
# Constants
rpP = 0.0625
rpFe = 6.25e-5
k = 0.1
ref_PN = 1.04
ref_FeN = 1.2

#bio_PN = (P_tot/N_tot)*(1/rpP)
#bio_FeN = (Fe_tot/N_tot)*(1/rpFe)

bio_PN = (P_remin/N_remin)*(1/rpP)
bio_FeN = (Fe_other/N_remin)*(1/rpFe)

#%% Make seasonal P:N and Fe:N arrays

season = ['DJF','MAM','JJA','SON']
DJF = [range(-1,1),range(2,4),range(5,7),range(8,10)]

diaz_int_DJF = (diaz_int[0,:,:]+diaz_int[1,:,:]+diaz_int[11,:,:])/3
bio_PN_DJF = (bio_PN[0,:,:]+bio_PN[1,:,:]+bio_PN[11,:,:])/3
bio_FeN_DJF = (bio_FeN[0,:,:]+bio_FeN[1,:,:]+bio_FeN[11,:,:])/3

#%% mask where P:N OR Fe:N is sufficient to support diazotrophs
mask_JJA = np.where((np.mean(bio_FeN[DJF[2],:,:],axis=0) > ref_FeN) & (np.mean(bio_PN[DJF[2],:,:],axis=0) > ref_PN), 1, 0)

#%% Plot seasonal nutrient ratio patterns

levs_phi = np.linspace(-1.2,-0.8,11)

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(4,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,12),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,4):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.05,0.9,''+(str(season[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    #con = ax[i].contour(lon,lat,np.mean(bio_PN[:,:,:],axis=0),levels=[1.04],colors='r',linewidths=1,linstyle='solid')
c1 = ax[0].contourf(lon,lat,bio_PN_DJF[:,:],levels=levs_phi,cmap=colmap,extend='both')
c2 = ax[1].contourf(lon,lat,np.mean(bio_PN[2:4,:,:],axis=0),levels=levs_phi,cmap=colmap,extend='both')
c3 = ax[2].contourf(lon,lat,np.mean(bio_PN[5:7,:,:],axis=0),levels=levs_phi,cmap=colmap,extend='both')
c4 = ax[3].contourf(lon,lat,np.mean(bio_PN[8:10,:,:],axis=0),levels=levs_phi,cmap=colmap,extend='both')
cbar1 = plt.colorbar(c1,ax=ax[0])
cbar2 = plt.colorbar(c2,ax=ax[1])
cbar3 = plt.colorbar(c3,ax=ax[2])
cbar4 = plt.colorbar(c4,ax=ax[3])
ax[0].set_title('Seasonal nutrient ratio P:N (remin only)')
#cbar.set_label(''+str(name_nut[nu])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/seasonal_PN_remin.png', bbox_inches='tight', dpi=300)

#%%
levs_phi = np.linspace(-2,1,11)

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(4,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,12),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,4):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.05,0.9,''+(str(season[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    #con = ax[i].contour(lon,lat,np.mean(bio_FeN[:,:,:],axis=0),levels=[1.2],colors='r',linewidths=1,linstyle='solid')
c1 = ax[0].contourf(lon,lat,bio_FeN_DJF[:,:],levels=levs_phi,cmap=colmap,extend='both')
c2 = ax[1].contourf(lon,lat,np.mean(bio_FeN[2:4,:,:],axis=0),levels=levs_phi,cmap=colmap,extend='both')
c3 = ax[2].contourf(lon,lat,np.mean(bio_FeN[5:7,:,:],axis=0),levels=levs_phi,cmap=colmap,extend='both')
c4 = ax[3].contourf(lon,lat,np.mean(bio_FeN[8:10,:,:],axis=0),levels=levs_phi,cmap=colmap,extend='both')
cbar1 = plt.colorbar(c1,ax=ax[0])
cbar2 = plt.colorbar(c2,ax=ax[1])
cbar3 = plt.colorbar(c3,ax=ax[2])
cbar4 = plt.colorbar(c4,ax=ax[3])
ax[0].set_title('Seasonal nutrient ratio Fe:N')
#cbar.set_label(''+str(name_nut[nu])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/seasonal_FeN_remin.png', bbox_inches='tight', dpi=300)

#%% Plot 1 nutrient at 1 depth

levs_diaz = np.linspace(0,20,21)

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(4,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,12),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,4):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.05,0.9,''+(str(season[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c1 = ax[0].contourf(lon,lat,diaz_int_DJF[:,:],levels=levs_diaz,cmap=colmap,extend='both')
c2 = ax[1].contourf(lon,lat,np.mean(diaz_int[2:4,:,:],axis=0),levels=levs_diaz,cmap=colmap,extend='both')
c3 = ax[2].contourf(lon,lat,np.mean(diaz_int[5:7,:,:],axis=0),levels=levs_diaz,cmap=colmap,extend='both')
c4 = ax[3].contourf(lon,lat,np.mean(diaz_int[8:10,:,:],axis=0),levels=levs_diaz,cmap=colmap,extend='both')
cbar1 = plt.colorbar(c1,ax=ax[0])
cbar2 = plt.colorbar(c2,ax=ax[1])
cbar3 = plt.colorbar(c3,ax=ax[2])
cbar4 = plt.colorbar(c4,ax=ax[3])
ax[0].set_title('Seasonal diazotroph biogeography')
#cbar.set_label(''+str(name_nut[nu])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/seasonal_diaz.png', bbox_inches='tight', dpi=300)


#%% Plot 1 nutrient at 1 depth
colmap = plt.get_cmap('RdBu_r')
levs = np.linspace(-10,10,21)

diaz_mean = np.mean(diaz_int[:,:,:],axis=0)

fig,ax = plt.subplots(4,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)}, sharex=True, sharey=True, figsize=(12,7))
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
m = 0
for i in range(4):
    for j in range(3):
        #ax[i,j].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
        ax[i,j].coastlines(color='#888888',linewidth=1.5)
        ax[i,j].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
        c0 = ax[i,j].contourf(lon,lat,(diaz_mean[:,:]-diaz_int[m,:,:]),levels=levs,extend='both',cmap=colmap)
        #c0 = ax[i,j].contourf(lon,lat,diaz_int[m,:,:],levels=levs,extend='both',cmap=colmap)
        #ax[i,j].set_extent([220-360, 295-360, -30, 30])
        #ax[i,j].set_xticks([240, 270], crs=ccrs.PlateCarree())
        #ax[i,j].set_yticks([-20, 0, 20], crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax[2,j].xaxis.set_major_formatter(lon_formatter)
        ax[i,0].yaxis.set_major_formatter(lat_formatter)       
        ax[i,j].grid()
        ax[i,j].text(0.78,0.9, months[m], size=10, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'),transform=ax[i,j].transAxes)
        m += 1    
        plt.suptitle('Monthly anomalies from mean diazotroph biogeography',y=0.92)
        fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
        cbar = fig.colorbar(c0, cax=cbar_ax)
        #cbar.set_label('more         (m)           less', rotation=90,labelpad=10)
#fig.savefig('/Users/meilers/MITinternship/Plots/monthly_diaz_anomalies2.png', bbox_inches='tight', dpi=300)
        


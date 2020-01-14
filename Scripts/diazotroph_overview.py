# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cmocean as cm
import matplotlib.cm as cmo
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from copy import deepcopy
import scipy.interpolate as si

#%%############################################################################
############### Load data from Darwin model output ############################
###############################################################################

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/supply50m.nc')
area_info = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')

lon = grid.lon   #needed for plotting
lat = grid.lat    #needed for plotting

area = area_info.rA

# Read in diazotroph data from model - same setup as for nutrients
# Simulation: run19_33
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

#%%############################################################################
################# MAREDAT - Load and condense diazotroph data #################
###############################################################################

ds = xr.open_dataset('/Users/meilers/MITinternship/Data/MarEDat20130403Diazotrophs.nc',decode_times=False)

# extract variables which are needed and convert/integrate
lon_d = ds['LONGITUDE']
lat_d = ds['LATITUDE'][10:-10] # to match the latitude of the nutrient data

obs = ds['OBSERVATIONS']
abund = ds['ABUNDANCE']
bm = ds['BIOMASS']
nifH = ds['nifHbiom']

obs_tot = np.sum(obs[:,:,10:-10,:],axis=(0,1))
abund_tot = np.sum(abund[:,:,10:-10,:],axis=(0,1))
bm_tot = np.sum(bm[:,:,:,:],axis=(0,1))
nifH_tot = np.sum(nifH[:,:,10:-10,:],axis=(0,1))

#%%############################################################################
################### Tang and Cassar database ##################################
###############################################################################

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
#longitude = pd.DataFrame(diazotroph_observations, columns = ['LONGITUDE'])
#latitude = pd.DataFrame(diazotroph_observations, columns = ['LATITUDE'])

#nifH_Tri = pd.DataFrame(diazotroph_observations, columns = ['Trichodesmium nifH Gene (x106 copies m-2)'])
#nifH_UCYN_A = pd.DataFrame(diazotroph_observations, columns = ['UCYN-A nifH Gene (x106 copies m-2)'])
#nifH_UCYN_B = pd.DataFrame(diazotroph_observations, columns = ['UCYN-B nifH Gene (x106 copies m-2)'])
#nifH_UCYN_C = pd.DataFrame(diazotroph_observations, columns = ['UCYN-C nifH Gene (x106 copies m-2)'])
#nifH_Richelia = pd.DataFrame(diazotroph_observations, columns = ['Richelia nifH Gene (x106 copies m-2)'])
#nifH_Calothrix = pd.DataFrame(diazotroph_observations, columns = ['Calothrix nifH Gene  (x106 copies m-2)'])
#nifH_Gamma = pd.DataFrame(diazotroph_observations, columns = ['Gamma nifH Gene (x106 copies/m3)'])

# Single columns of nifH database in list format
#nifH_Tri = diazotroph_observations['Trichodesmium nifH Gene (x106 copies m-2)'].values.tolist()
nifH_Tri = diazotroph_observations['Trichodesmium nifH Gene (x106 copies m-2)'].astype(np.float32)
nifH_UCYN_A = diazotroph_observations['UCYN-A nifH Gene (x106 copies m-2)'].astype(np.float32)
nifH_UCYN_B = diazotroph_observations['UCYN-B nifH Gene (x106 copies m-2)'].astype(np.float32)
nifH_UCYN_C = diazotroph_observations['UCYN-C nifH Gene (x106 copies m-2)'].astype(np.float32)
nifH_Richelia = diazotroph_observations['Richelia nifH Gene (x106 copies m-2)'].astype(np.float32)
nifH_Calothrix = diazotroph_observations['Calothrix nifH Gene  (x106 copies m-2)'].astype(np.float32)
nifH_Gamma = diazotroph_observations['Gamma nifH Gene (x106 copies/m3)'].astype(np.float32)
lon_nifH = diazotroph_observations['LONGITUDE'].astype(np.float32)
lat_nifH = diazotroph_observations['LATITUDE'].astype(np.float32)
year = diazotroph_observations['YEAR'].astype(np.float32)
month = diazotroph_observations['MONTH'].astype(np.float32)
#nifH = nifH_Tri + nifH_UCYN_A + nifH_UCYN_B + nifH_UCYN_C + nifH_Richelia + nifH_Calothrix + nifH_Gamma
#nifH_sum = np.sum(nifH, axis=(1))

#%%compile the nifH data of different species into 1 matrix
var_list = [lon_nifH, lat_nifH, year, month, nifH_Tri, nifH_UCYN_A, nifH_UCYN_B, nifH_UCYN_C, nifH_Richelia, nifH_Calothrix, nifH_Gamma]
nifH_matrix = np.zeros((len(nifH_Tri),len(var_list)+1))
for i in range(len(var_list)):
    nifH_matrix[:,i] = var_list[i]

# find datapoints where no nifH gene is found at all (loop over all species)
for j in range(len(nifH_Tri)):
    if np.nansum(nifH_matrix[j,4:11]) > 0:
        nifH_matrix[j,-1] = 1

# create list of indices of nifH presences and absences. Can be used for plotting or further data manipulation
presence = np.where(nifH_matrix[:,-1] > 0)
absence = np.where(nifH_matrix[:,-1] == 0)

# get lon, lat of presence and absence
lon_pres = lon_nifH[presence[0]].astype(np.float32)
lat_pres = lat_nifH[presence[0]].astype(np.float32)
lon_abs = lon_nifH[absence[0]]
lat_abs = lat_nifH[absence[0]]

#%% create inteprolator
grid = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')
lon = grid.X.values
lat = grid.Y.values
vals = grid.XC.values  # for testing - replace by data to be interpolated

I = si.RegularGridInterpolator((lat, lon), vals, 'nearest',
                               bounds_error=False, fill_value=None)

# where to evaluate data
#lon_d = np.r_[0., 10., 20.]
#lat_d = np.r_[-80., 0., 10.]

n_d = len(lon_pres)
latlon = np.zeros((n_d, 2))
latlon[:,0] = lat_pres
latlon[:,1] = np.mod(lon_pres, 360.)

# interpolate
vals_d = I(latlon)

#%% Presence/absence on monthly time scales
# find a way to display the data on monthly scales
pres_month = nifH_matrix[presence,3]
abs_month = nifH_matrix[absence,3]


#%% Import lat, lon, and regions file to define ocean basins/regions
# after Teng et al., 2014

regions = pd.read_csv(r'/Users/meilers/MITinternship/Data/Regions/regions.csv')
reg_lon = pd.read_csv(r'/Users/meilers/MITinternship/Data/Regions/lons.csv')
reg_lat = pd.read_csv(r'/Users/meilers/MITinternship/Data/Regions/lats.csv')


#%% Manipulating the nifH data to bring it into mappable form - SPECIES SPECIFIC

Tri_list = np.where(nifH_Tri > 0)
UCYN_A_list = np.where(nifH_UCYN_A > 0) 
UCYN_B_list = np.where(nifH_UCYN_B > 0)
UCYN_C_list = np.where(nifH_UCYN_C > 0)
Richelia_list = np.where(nifH_Richelia > 0)
Calothrix_list = np.where(nifH_Calothrix > 0)
Gamma_list = np.where(nifH_Gamma > 0)
#lon_nifH[Tri_list[0]]
#lat_nifH[Tri_list[0]]

#%% Absence of nifH (zeros or n.d.)

no_Tri_list = np.where(nifH_Tri == 0)
no_UCYN_A_list = np.where(nifH_UCYN_A == 0) 
no_UCYN_B_list = np.where(nifH_UCYN_B == 0)
no_UCYN_C_list = np.where(nifH_UCYN_C == 0)
no_Richelia_list = np.where(nifH_Richelia == 0)
no_Calothrix_list = np.where(nifH_Calothrix == 0)
no_Gamma_list = np.where(nifH_Gamma == 0)

#%% define constants and dz
sec = 1  #seconds per year (365.25*86400)
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
    #print(np.max(diaz_int[:,i,:,:]))
    #print(i)
diaz_int = np.sum(diaz_int,axis=1)
diaz_std = deepcopy(diaz_int)
diaz_std = np.std(diaz_std,axis=0)
diaz_int = np.mean(diaz_int,axis=0)
diaz_cv = diaz_std/diaz_int

#diaz_int_100 = np.zeros((12,6,160,360))
#for i in range(len(dz)):
#    diaz_int_100[:,i,:,:] = diaz[:,i,:,:]*dz[i]*sec  
#    #print(np.max(diaz_int_100[:,i,:,:]))
#    #print(i)
#diaz_int_100 = np.sum(diaz_int_100,axis=1)
#diaz_int_100 = np.mean(diaz_int_100,axis=0)

#%% mask where diazotroph biomass is simulated
mask = np.where((diaz_int > 1e-04), 1, 0)
mask_out = np.where((diaz_int < 1e-04), 1, 0)

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

#%% Calculate accuracy for absences
list_idx = 4
lat_corr_abs = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
lon_corr_abs = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
OUT = np.sum(mask_out[lat_corr_abs,lon_corr_abs])/len(lat_corr_abs)
print(OUT)


#%% Same IN/OUT calculation for Tang and Cassar data

#lat_corr = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
#lon_corr = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
IN_nifH = np.sum(mask[lat_nifH,lon_nifH])/len(lat_nifH)
print(IN_nifH)

#%% Calculate accuracy for absences
list_idx = 4
lat_corr_abs = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
lon_corr_abs = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
OUT = np.sum(mask_out[lat_corr_abs,lon_corr_abs])/len(lat_corr_abs)
print(OUT)


#%% Plot diazotroph biomass simulated in Darwin - 1 subplot (mean, STD, or CV)

col = cm.cm.haline

which_index = 2 #chose which metric to plot: 0 = mean, 1 = std, 2 = cv
which_metric = [diaz_int, diaz_std, diaz_cv]
which_level = [np.linspace(0,40,21),np.linspace(0,20,21),np.linspace(0,5,21)]
which_text = ['mean','std','CV']
which_label = ['mmolC m$^{-2}$','mmolC m$^{-2}$','[-]']


fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
ax.text(0.2,0.9,''+str(which_text[which_index])+'',transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,which_metric[which_index],levels=which_level[which_index],cmap=col,extend='max')
#ax.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='orange',label='presence')
#ax.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='m',label='absence')
#ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.015, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label(''+str(which_label[which_index])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_'+str(which_text[which_index])+'.png', bbox_inches='tight', dpi=300)

#%% Plot diazotroph biomass simulated in Darwin (integrated over entire depth range; averaged over 1 year; STD; coefficient variation CV)

#col = plt.get_cmap('RdBu_r')
col = cm.cm.haline

depth_lab = ['mean','std','CV']
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
    ax[i].text(0.2,0.9,''+(str(depth_lab[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax[0].contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
c1 = ax[1].contourf(lon,lat,diaz_std,levels=np.linspace(0,20,21),cmap=col,extend='max')
c2 = ax[2].contourf(lon,lat,diaz_cv,levels=np.linspace(0,5,21),cmap=col,extend='max')
cbar0 = plt.colorbar(c0,ax=ax[0])
cbar1 = plt.colorbar(c1,ax=ax[1])
cbar2 = plt.colorbar(c2,ax=ax[2])
cbar0.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
cbar1.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
cbar2.set_label('[-]',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nodata_alldepth.png', bbox_inches='tight', dpi=300)

#%% Plot diazotroph biomass simulated in Darwin and Tang data for nifH gene counts

col = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.text(0.2,0.9,''+(str(depth_lab[1])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
#ax.plot(lon_nifH[Tri_list[0]],lat_nifH[Tri_list[0]],'.',color='orange',label='Trichodesmium')
#ax.plot(lon_nifH[UCYN_A_list[0]],lat_nifH[UCYN_A_list[0]],'.',color='orange',label='UCYN-A')
#ax.plot(lon_nifH[UCYN_B_list[0]],lat_nifH[UCYN_B_list[0]],'.',color='orange',label='UCYN-B')
#ax.plot(lon_nifH[UCYN_C_list[0]],lat_nifH[UCYN_C_list[0]],'.',color='orange',label='UCYN-C')
#ax.plot(lon_nifH[Richelia_list[0]],lat_nifH[Richelia_list[0]],'.',color='orange',label='Richelia')
#ax.plot(lon_nifH[Calothrix_list[0]],lat_nifH[Calothrix_list[0]],'.',color='orange',label='Calothrix')
ax.plot(lon_nifH[Gamma_list[0]],lat_nifH[Gamma_list[0]],'.',color='orange',label='Gamma')
#ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifHpresence.png', bbox_inches='tight', dpi=300)

#%% Plot diazotroph biomass simulated in Darwin and Tang data for nifH ABSENCES

col = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.text(0.2,0.9,''+(str(depth_lab[1])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
ax.plot(lon_nifH[no_Tri_list[0]],lat_nifH[no_Tri_list[0]],'x',color='m',label='Trichodesmium')
ax.plot(lon_nifH[no_UCYN_A_list[0]],lat_nifH[no_UCYN_A_list[0]],'x',color='m',label='UCYN-A')
ax.plot(lon_nifH[no_UCYN_B_list[0]],lat_nifH[no_UCYN_B_list[0]],'x',color='m',label='UCYN-B')
ax.plot(lon_nifH[no_UCYN_C_list[0]],lat_nifH[no_UCYN_C_list[0]],'x',color='m',label='UCYN-C')
ax.plot(lon_nifH[no_Richelia_list[0]],lat_nifH[no_Richelia_list[0]],'x',color='m',label='Richelia')
ax.plot(lon_nifH[no_Calothrix_list[0]],lat_nifH[no_Calothrix_list[0]],'x',color='m',label='Calothrix')
ax.plot(lon_nifH[no_Gamma_list[0]],lat_nifH[no_Gamma_list[0]],'x',color='m',label='Gamma')
#ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifHabsence.png', bbox_inches='tight', dpi=300)

#%%############################################################################ 
################### Best figure for now #######################################
###############################################################################

#col = plt.get_cmap('RdBu_r')
col = cm.cm.haline

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.text(0.2,0.9,''+(str(depth_lab[1])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
ax.plot(lon_nifH[presence[0]],lat_nifH[presence[0]],'.',color='orange',label='any nifH present')#label='Trichodesmium')
ax.plot(lon_nifH[absence[0]],lat_nifH[absence[0]],'x',color='m',label='all nifH absent')
ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifH.png', bbox_inches='tight', dpi=300)

#%% 
# side-by-side plot with one map only having absences (maybe colour by month 
# of observation?), and the other having only presence of nifH

#col = plt.get_cmap('RdBu_r')
#col = cm.cm.haline
#
#fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
#lon_formatter = LongitudeFormatter(zero_direction_label=True)
#lat_formatter = LatitudeFormatter()
#ax.coastlines(color='#888888',linewidth=1.5)
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)
#ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
#ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
##ax.text(0.2,0.9,''+(str(depth_lab[1])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#c0 = ax.contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
#ax.plot(lon_nifH[presence[0]],lat_nifH[presence[0]],'.',color='orange',label='any nifH present')#label='Trichodesmium')
#ax.plot(lon_nifH[absence[0]],lat_nifH[absence[0]],'x',color='m',label='all nifH absent')
#ax.legend(loc='best')
#fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
#cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
#cbar = fig.colorbar(c0, cax=cbar_ax)
#cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
#plt.show()
jan_pres = np.where(nifH_matrix[presence][:,3] == 1)
feb_pres = np.where(nifH_matrix[presence][:,3] == 2)
mar_pres = np.where(nifH_matrix[presence][:,3] == 3)
apr_pres = np.where(nifH_matrix[presence][:,3] == 4)
may_pres = np.where(nifH_matrix[presence][:,3] == 5)
jun_pres = np.where(nifH_matrix[presence][:,3] == 6)
jul_pres = np.where(nifH_matrix[presence][:,3] == 7)
aug_pres = np.where(nifH_matrix[presence][:,3] == 8)
sep_pres = np.where(nifH_matrix[presence][:,3] == 9)
oct_pres = np.where(nifH_matrix[presence][:,3] == 10)
nov_pres = np.where(nifH_matrix[presence][:,3] == 11)
dec_pres = np.where(nifH_matrix[presence][:,3] == 12)

jan_abs = np.where(nifH_matrix[absence][:,3] == 1)
feb_abs = np.where(nifH_matrix[absence][:,3] == 2)
mar_abs = np.where(nifH_matrix[absence][:,3] == 3)
apr_abs = np.where(nifH_matrix[absence][:,3] == 4)
may_abs = np.where(nifH_matrix[absence][:,3] == 5)
jun_abs = np.where(nifH_matrix[absence][:,3] == 6)
jul_abs = np.where(nifH_matrix[absence][:,3] == 7)
aug_abs = np.where(nifH_matrix[absence][:,3] == 8)
sep_abs = np.where(nifH_matrix[absence][:,3] == 9)
oct_abs = np.where(nifH_matrix[absence][:,3] == 10)
nov_abs = np.where(nifH_matrix[absence][:,3] == 11)
dec_abs = np.where(nifH_matrix[absence][:,3] == 12)

#%%
colors = cmo.rainbow(np.linspace(0, 1, 12))
depth_lab = ['presences','absences']
howmany = [len(presence[0]),len(absence[0])]
fig,ax = plt.subplots(2,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,6),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,2):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.75,0.9,''+(str(depth_lab[i])+' n='+str(howmany[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax[0].contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
c1 = ax[1].contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
ax[0].plot(nifH_matrix[presence][:,0][jan_pres][0],nifH_matrix[presence][:,1][jan_pres][0],'.',color=colors[0],label='Jan')
ax[0].plot(nifH_matrix[presence][:,0][feb_pres[0]],nifH_matrix[presence][:,1][feb_pres[0]],'.',color=colors[1],label='Feb')
ax[0].plot(nifH_matrix[presence][:,0][mar_pres[0]],nifH_matrix[presence][:,1][mar_pres[0]],'.',color=colors[2],label='Mar')
ax[0].plot(nifH_matrix[presence][:,0][apr_pres[0]],nifH_matrix[presence][:,1][apr_pres[0]],'.',color=colors[3],label='Apr')
ax[0].plot(nifH_matrix[presence][:,0][may_pres[0]],nifH_matrix[presence][:,1][may_pres[0]],'.',color=colors[4],label='May')
ax[0].plot(nifH_matrix[presence][:,0][jun_pres[0]],nifH_matrix[presence][:,1][jun_pres[0]],'.',color=colors[5],label='Jun')
ax[0].plot(nifH_matrix[presence][:,0][jul_pres[0]],nifH_matrix[presence][:,1][jul_pres[0]],'.',color=colors[6],label='Jul')
ax[0].plot(nifH_matrix[presence][:,0][aug_pres[0]],nifH_matrix[presence][:,1][aug_pres[0]],'.',color=colors[7],label='Aug')
ax[0].plot(nifH_matrix[presence][:,0][sep_pres[0]],nifH_matrix[presence][:,1][sep_pres[0]],'.',color=colors[8],label='Sep')
ax[0].plot(nifH_matrix[presence][:,0][oct_pres[0]],nifH_matrix[presence][:,1][oct_pres[0]],'.',color=colors[9],label='Oct')
ax[0].plot(nifH_matrix[presence][:,0][nov_pres[0]],nifH_matrix[presence][:,1][nov_pres[0]],'.',color=colors[10],label='Nov')
ax[0].plot(nifH_matrix[presence][:,0][dec_pres[0]],nifH_matrix[presence][:,1][dec_pres[0]],'.',color=colors[11],label='Dec')
ax[1].plot(nifH_matrix[absence][:,0][jan_abs[0]],nifH_matrix[absence][:,1][jan_abs[0]],'x',color=colors[0],label='Jan')
ax[1].plot(nifH_matrix[absence][:,0][feb_abs[0]],nifH_matrix[absence][:,1][feb_abs[0]],'x',color=colors[1],label='Feb')
ax[1].plot(nifH_matrix[absence][:,0][mar_abs[0]],nifH_matrix[absence][:,1][mar_abs[0]],'x',color=colors[2],label='Mar')
ax[1].plot(nifH_matrix[absence][:,0][apr_abs[0]],nifH_matrix[absence][:,1][apr_abs[0]],'x',color=colors[3],label='Apr')
ax[1].plot(nifH_matrix[absence][:,0][may_abs[0]],nifH_matrix[absence][:,1][may_abs[0]],'x',color=colors[4],label='May')
ax[1].plot(nifH_matrix[absence][:,0][jun_abs[0]],nifH_matrix[absence][:,1][jun_abs[0]],'x',color=colors[5],label='Jun')
ax[1].plot(nifH_matrix[absence][:,0][jul_abs[0]],nifH_matrix[absence][:,1][jul_abs[0]],'x',color=colors[6],label='Jul')
ax[1].plot(nifH_matrix[absence][:,0][aug_abs[0]],nifH_matrix[absence][:,1][aug_abs[0]],'x',color=colors[7],label='Aug')
ax[1].plot(nifH_matrix[absence][:,0][sep_abs[0]],nifH_matrix[absence][:,1][sep_abs[0]],'x',color=colors[8],label='Sep')
ax[1].plot(nifH_matrix[absence][:,0][oct_abs[0]],nifH_matrix[absence][:,1][oct_abs[0]],'x',color=colors[9],label='Oct')
ax[1].plot(nifH_matrix[absence][:,0][nov_abs[0]],nifH_matrix[absence][:,1][nov_abs[0]],'x',color=colors[10],label='Nov')
ax[1].plot(nifH_matrix[absence][:,0][dec_abs[0]],nifH_matrix[absence][:,1][dec_abs[0]],'x',color=colors[11],label='Dec')
#ax[1].legend(loc='best')
cbar0 = plt.colorbar(c0,ax=ax[0])
cbar1 = plt.colorbar(c1,ax=ax[1])
cbar0.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
cbar1.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.legend(ncol=6,loc='lower center',bbox_to_anchor=(0.5,-0.4))
plt.show()
fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifH_pres-abs2.png', bbox_inches='tight', dpi=300)

#%% Plot diazotroph biomass simulated in Darwin and Tang data for nifH nifH gene counts & ABSENCES

#col = plt.get_cmap('RdBu_r')
col = cm.cm.haline

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.text(0.2,0.9,''+(str(depth_lab[1])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
c0 = ax.contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=col,extend='max')
ax.plot(lon_nifH[no_Tri_list[0]],lat_nifH[no_Tri_list[0]],'x',color='m',label='nifH absence')#label='Trichodesmium')
ax.plot(lon_nifH[no_UCYN_A_list[0]],lat_nifH[no_UCYN_A_list[0]],'x',color='m')#,label='UCYN-A')
ax.plot(lon_nifH[no_UCYN_B_list[0]],lat_nifH[no_UCYN_B_list[0]],'x',color='m')#,label='UCYN-B')
ax.plot(lon_nifH[no_UCYN_C_list[0]],lat_nifH[no_UCYN_C_list[0]],'x',color='m')#,label='UCYN-C')
ax.plot(lon_nifH[no_Richelia_list[0]],lat_nifH[no_Richelia_list[0]],'x',color='m')#,label='Richelia')
ax.plot(lon_nifH[no_Calothrix_list[0]],lat_nifH[no_Calothrix_list[0]],'x',color='m')#,label='Calothrix')
ax.plot(lon_nifH[no_Gamma_list[0]],lat_nifH[no_Gamma_list[0]],'x',color='m')#,label='Gamma')
ax.plot(lon_nifH[Tri_list[0]],lat_nifH[Tri_list[0]],'.',color='orange',label='nifH presence')#label='Trichodesmium')
ax.plot(lon_nifH[UCYN_A_list[0]],lat_nifH[UCYN_A_list[0]],'.',color='orange')#,label='UCYN-A')
ax.plot(lon_nifH[UCYN_B_list[0]],lat_nifH[UCYN_B_list[0]],'.',color='orange')#,label='UCYN-B')
ax.plot(lon_nifH[UCYN_C_list[0]],lat_nifH[UCYN_C_list[0]],'.',color='orange')#,label='UCYN-C')
ax.plot(lon_nifH[Richelia_list[0]],lat_nifH[Richelia_list[0]],'.',color='orange')#,label='Richelia')
ax.plot(lon_nifH[Calothrix_list[0]],lat_nifH[Calothrix_list[0]],'.',color='orange')#,label='Calothrix')
ax.plot(lon_nifH[Gamma_list[0]],lat_nifH[Gamma_list[0]],'.',color='orange')#,label='Gamma')
ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifH.png', bbox_inches='tight', dpi=300)

#%% Mask from Darwin biomass and nifH presence/absence
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(12,4))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,mask,levels=np.linspace(0,1,11),cmap=cm.cm.haline,extend='both')
#ax.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='orange',label='presence')
#ax.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='m',label='absence')
ax.plot(lon_nifH[Tri_list[0]],lat_nifH[Tri_list[0]],'.',color='orange',label='presence')
ax.plot(lon_nifH[UCYN_A_list[0]],lat_nifH[UCYN_A_list[0]],'.',color='orange')
ax.plot(lon_nifH[UCYN_B_list[0]],lat_nifH[UCYN_B_list[0]],'.',color='orange')
ax.plot(lon_nifH[UCYN_C_list[0]],lat_nifH[UCYN_C_list[0]],'.',color='orange')
ax.plot(lon_nifH[Richelia_list[0]],lat_nifH[Richelia_list[0]],'.',color='orange')
ax.plot(lon_nifH[Calothrix_list[0]],lat_nifH[Calothrix_list[0]],'.',color='orange')
ax.plot(lon_nifH[Gamma_list[0]],lat_nifH[Gamma_list[0]],'.',color='orange')
ax.plot(lon_nifH[no_Tri_list[0]],lat_nifH[no_Tri_list[0]],'.',color='m',label='absence')
ax.plot(lon_nifH[no_UCYN_A_list[0]],lat_nifH[no_UCYN_A_list[0]],'.',color='m')
ax.plot(lon_nifH[no_UCYN_B_list[0]],lat_nifH[no_UCYN_B_list[0]],'.',color='m')
ax.plot(lon_nifH[no_UCYN_C_list[0]],lat_nifH[no_UCYN_C_list[0]],'.',color='m')
ax.plot(lon_nifH[no_Richelia_list[0]],lat_nifH[no_Richelia_list[0]],'.',color='m')
ax.plot(lon_nifH[no_Calothrix_list[0]],lat_nifH[no_Calothrix_list[0]],'.',color='m')
ax.plot(lon_nifH[no_Gamma_list[0]],lat_nifH[no_Gamma_list[0]],'.',color='m')
ax.legend(loc='best')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#cbar = plt.colorbar(c,ax=ax)
#cbar.set_label(''+str(name_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_mask_nifHpres-abs.png', bbox_inches='tight', dpi=300)


#%% Plot diazotroph biomass from MAREDAT

col = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,4))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-75, -60, -40, -20, 0, 20, 40, 60, 75], crs=ccrs.PlateCarree())
ax.text(0.15,0.9,''+(str(depth_lab[0])+''),transform=ax.transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#c0 = ax.contourf(lon,lat,bm_tot,levels=np.linspace(0,4,21),cmap=cm.cm.haline)#,extend='both')
ax.scatter(lon_d[bm_tot[1]],lat_d[bm_tot[0]],'.',color='orange',label='presence')
#ax.plot(lon_d[absent_obs[1]],lat_d[absent_obs[0]],'.',color='m',label='absence')
ax.legend(loc='best')
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_withdata.png', bbox_inches='tight', dpi=300)

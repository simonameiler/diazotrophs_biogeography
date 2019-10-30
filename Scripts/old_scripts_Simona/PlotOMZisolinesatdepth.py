#!/usr/bin/env python

'''
Compute Hovmoller data for tracer and offshore transport from monthly mean output
and eventually store it in a pickle file for plotting.

Martin Frischknecht, Nov 2016
Adapted by Simona Meiler, June 2019
'''

#########################################################
# Load requirements and set switches
#########################################################

# Import python modules
import numpy as np
import xarray as xr
import netCDF4
import matplotlib.pyplot as plt
import scipy.stats as sp
from scipy.interpolate import interp1d
import cmocean as cm
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
exec(open('/home/koehne/Documents/scripts/python/martinfr/myromstools.py').read()) # import compute_zlev, compute_dz

#%% Read grid
gridfile = '/net/kryo/work/koehne/roms/inputs/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/grd/humpac15_grd.nc'
romsGrd = getGrid(gridfile)
# Add necessary attributes
filepath = '/net/kryo/work/koehne/roms/output/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_hindcast_1979_2016_6/monthly/his/'
romsGrd.getAttrs(filepath+'humpac15_1979_his.nc')
# Add lat lon
romsGrd.getLatLon()
# Add grid area
romsGrd.getArea()
# Add angle
romsGrd.getAngle()

grid = xr.open_dataset(gridfile)

#%% Run
setup = 'humpac15'
runtag = '{}_hindcast_1979_2016_6'.format(setup)
inpath = '/net/kryo/work/martinfr/Roms/Output/humpac15/{}/avg/'.format(runtag)
inpath = '/net/kryo/work/koehne/roms/output/{}/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/{}//monthly/avg/'.format(setup,runtag)
# Set path to save output
outpath = '/net/kryo/work/meilers/Data/{}/{}/'.format(setup,runtag)
# Filelist
year_start = 1979
year_end = 2017
filelist = [inpath+'humpac15_{}_avg.nc'.format(yr) for yr in range(year_start,year_end)]
outfile = outpath+'OMZ_depths_1979_2017_interpolated.nc'

#%%######################################################
# Plot
#########################################################
nc_out = netCDF4.Dataset(outfile,'r')
#% what to plot

# read in Niño 3.4 index from model (not observed)
nino34_idx = np.load('/home/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/ENSOindices/Nino_indices/Nino34_model.npy')
# read in the normalized and smoothed (5-month running mean) Niño 3.4 index (from model)
nino34_idx_norm = np.load('/home/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/ENSOindices/Nino_indices/Nino34_model_norm.npy')
# define El Niño months, where anomalies > 0.4 (source Martin Frischknecht)

m = 3
kernel = np.ones(m)/m #to calculate the ONI we use a 3-months running mean (over the modeled Niño 3.4 anomalies; meaning: SST-30-year climatology)
oni_anom = np.convolve(nino34_idx,kernel,mode='same')

# this block is only used for plotting
from copy import deepcopy
oni_anom_pos = deepcopy(oni_anom)
oni_anom_neg = deepcopy(oni_anom)
oni_anom_pos[oni_anom_pos<0.5]=0
oni_anom_neg[oni_anom_neg>-0.5]=0

el_nino_months = list(np.where(oni_anom>0.5))
la_nina_months = list(np.where(oni_anom<-0.5))
#neutral_months = list(np.where(np.logical_and(oni_anom<=0.5, oni_anom>=-0.5)))

#test = len(el_nino_months[0]) + len(la_nina_months[0]) + len(neutral_months[0]) #to check whether we assigned all months with an ENSO status

data_to_plot = np.mean(nc_out['depthOMZ'][:,:,:,:],axis=0)
#neutral_to_plot= np.mean(nc_out['depthOMZ'][neutral_months[0],:,:,:], axis=0) # probably I won't use this
el_nino_to_plot = np.mean(nc_out['depthOMZ'][el_nino_months[0],:,:,:], axis=0)
la_nina_to_plot = np.mean(nc_out['depthOMZ'][la_nina_months[0],:,:,:], axis=0)

#%% Prepare data to plot

# Mask and thresholds
mask = np.zeros_like(grid.mask_rho) + np.NaN #creates a mask to fill all the continents with NaN
mask[grid.mask_rho==1] = 1

thresholds = [60,22,5]
pos = 0 #to chose which threshold to plot --> position in the list of thresholds
thresh = thresholds[pos] #to give out the value of the chosen threshold


#%% load masks
ETNP60mask3D = np.load('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/Masks/ETNP60mask3D.npy')
ETNP60mask = np.load('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/Masks/ETNP60mask2D.npy')
ETSP60mask3D = np.load('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/Masks/ETSP60mask3D.npy')
ETSP60mask = np.load('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/Masks/ETSP60mask2D.npy')
mask_100km3D = np.load('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/Masks/mask_100km3D.npy')
mask_100km = np.load('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/Masks/mask_100km.npy')

#%% Chose area over which to average
area_handle = ['ETNP','ETSP','ETSP coast 100km']
area_save = ['ETNP','ETSP','ETSP coast 100km']
area_pos = 0 #chose area here
masked_areas_3D = [ETNP60mask3D, ETSP60mask3D, mask_100km3D]
masked_areas_2D = [ETNP60mask, ETSP60mask, mask_100km]
mask_area_cho_3D = masked_areas_3D[area_pos]
mask_area_cho_2D = masked_areas_2D[area_pos]
#%% Data to plot

# Read data for different EN and LN events
EN97_to_plot = np.mean(nc_out['depthOMZ'][220:232,:,:,:], axis=0) #months 220:231 are El Niño months
EN82_to_plot = np.mean(nc_out['depthOMZ'][149:162,:,:,:], axis=0) #months 149:162 are El Niño months)
EN15_to_plot = np.mean(nc_out['depthOMZ'][436:447,:,:,:], axis=0) #months 436:447 are El Niño months)
LN88_to_plot = np.mean(nc_out['depthOMZ'][112:126,:,:,:], axis=0) #months 112:125 are La Niña months
EN02_to_plot = np.mean(nc_out['depthOMZ'][281:290,:,:,:], axis=0) #months 281:289 are El Niño months
LN07_to_plot = np.mean(nc_out['depthOMZ'][322:326,:,:,:], axis=0) #months 322:325 are La Niña months
LN95_to_plot = np.mean(nc_out['depthOMZ'][200:207,:,:,:], axis=0) #months 200:206 are La Niña months

diff_allyearsEN = data_to_plot[pos,:,:] - el_nino_to_plot[pos,:,:]
diff_allyearsLN = data_to_plot[pos,:,:] - la_nina_to_plot[pos,:,:]
diff_allyearsEN97 = data_to_plot[pos,:,:] - EN97_to_plot[pos,:,:] # EP El Niño event
diff_allyearsEN82 = data_to_plot[pos,:,:] - EN82_to_plot[pos,:,:] # EP El Niño event
diff_allyearsEN15 = data_to_plot[pos,:,:] - EN15_to_plot[pos,:,:] # EP El Niño event
diff_allyearsLN88 = data_to_plot[pos,:,:] - LN88_to_plot[pos,:,:] # EP La Niña event
diff_allyearsEN02 = data_to_plot[pos,:,:] - EN02_to_plot[pos,:,:] # CP El Niño event
diff_allyearsLN07 = data_to_plot[pos,:,:] - LN07_to_plot[pos,:,:] # CP La Niña event
diff_allyearsLN95 = data_to_plot[pos,:,:] - LN95_to_plot[pos,:,:] # CP La Niña event

# to show the relative change 
rel_changeEN = diff_allyearsEN/data_to_plot[pos,:,:]*100 #result in %
rel_changeLN = diff_allyearsLN/data_to_plot[pos,:,:]*100 #result in %
rel_changeEN97 = diff_allyearsEN97/data_to_plot[pos,:,:]*100 #result in %
rel_changeLN88 = diff_allyearsLN88/data_to_plot[pos,:,:]*100 #result in %
rel_changeEN02 = diff_allyearsEN02/data_to_plot[pos,:,:]*100 #result in %
rel_changeLN07 = diff_allyearsLN07/data_to_plot[pos,:,:]*100 #result in %
#test = rel_changeEN-rel_changeLN
#test_max = np.max(test)
#mean_of_them_all = np.mean(rel_changeEN)

#%% Plot the mean depth of the oxyisolines over all years and 3 thresholds
#maskeli = np.zeros_like(mask)
#mask[data_to_plot[0,:,:]<100] = False
##
fig,ax = plt.subplots()
#c0 = ax.contourf(grid.lon_rho,grid.lat_rho,mask)
ax.contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:]*mask,levels=[100],colors='w',linewidths=1)#,transform=ccrs.PlateCarree())
#ax.set_xlim([225,300])
#ax.set_ylim([-30,30])
ax.grid('on')

#%%
fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
ax[0].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
ax[1].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
ax[2].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
for i in range(len(thresholds)):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
    #ax[i].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
    c0 = ax[i].contourf(grid.lon_rho-180,grid.lat_rho,data_to_plot[i,:,:]*mask,levels=np.linspace(0,1000,11), cmap=cm.cm.ice, extend='max')
#    ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[i,:,:]*mask,levels=[100],colors='w',linewidths=1)
#    ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[i,:,:]*mask,levels=[400],linestyle='dashed',colors='w',linewidths=1)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[0].text(0.3,0.1,'60 mmol m$^{-3}$',transform=ax[0].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    ax[1].text(0.3,0.1,'22 mmol m$^{-3}$',transform=ax[1].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))          
    ax[2].text(0.3,0.1,'5 mmol m$^{-3}$',transform=ax[2].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))        
    ax[i].set_ylim([-60,80])
    ax[i].set_extent([180-360, 295-360, -40, 40])
#    ax[i].set_xticks([190, 220, 250, 280], crs=ccrs.PlateCarree())
#    ax[i].set_yticks([-40, -20, 0, 20, 40], crs=ccrs.PlateCarree())
#    ax[i].set_ylim([-40, 40])
    ax[i].grid()
#cb = plt.colorbar(c0,ax=ax[2], fraction=0.042, pad=0.04)
#ax[0].contour(grid.lon_rho-360.,grid.lat_rho,ETNP60mask,levels=[0.5],colors='w',linewidths=2,transform=ccrs.PlateCarree())
#ax[0].contour(grid.lon_rho-360.,grid.lat_rho,ETSP60mask,levels=[0.5],colors='w',linewidths=2,transform=ccrs.PlateCarree())
fig.subplots_adjust(right=0.85)
cb = fig.add_axes([0.9,0.25,0.01,0.5])
fig.colorbar(c0,cax=cb,label='Depth (m)')
#cb.set_label('[m]',rotation=90, position=(0.5,1))
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZisolinesatdepth_mean_'+str(year_start)+'-'+str(year_end)+'_allthresh.png', bbox_inches='tight', dpi=300)


#%% Plot the mean depth of the oxyisolines for different ENSO modes

EN = [diff_allyearsEN, diff_allyearsEN97, diff_allyearsEN02]
LN =[diff_allyearsLN, diff_allyearsLN88, diff_allyearsLN07]
ENnames = ['All El Niño','EN 97/98 (EP)','EN 02/03 (CP)']

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
for i in range(len(EN)):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
    c0 = ax[i].contourf(grid.lon_rho-180,grid.lat_rho,EN[i]*mask,levels=np.linspace(-100,100,11), cmap=colmap, extend='both')#,transform=ccrs.PlateCarree(central_longitude=180))
    c1 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[100],colors='k',linewidths=1)#,transform=ccrs.PlateCarree(central_longitude=180))
    c2 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[400],colors='k',linewidths=1,linestyles='dashed')#,transform=ccrs.PlateCarree(central_longitude=180))
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.3,0.1,''+str(ENnames[i])+'',transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    ax[i].set_ylim([-60,80])
    ax[i].set_extent([180-360, 295-360, -40, 40])
#    ax[i].set_xticks([190, 220, 250, 280], crs=ccrs.PlateCarree())
#    ax[i].set_yticks([-40, -20, 0, 20, 40], crs=ccrs.PlateCarree())
#    ax[i].set_ylim([-40, 40])
    ax[i].grid()
fig.subplots_adjust(right=0.85)
cb = fig.add_axes([0.9,0.25,0.01,0.5])
fig.colorbar(c0,cax=cb,label='m')
#cb.set_label('[m]',labelpad=-20, y=1.09, rotation=0)
#plt.subplots_adjust(wspace=0.05)
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZisolinesatdepth_mean_EN_'+str(thresh)+'.png', bbox_inches="tight", dpi=300)

#%% Plot to show the La Niña change of the depth of OMZ

EN = [diff_allyearsEN, diff_allyearsEN97, diff_allyearsEN02]
LN =[diff_allyearsLN, diff_allyearsLN95, diff_allyearsLN07]
ENnames = ['All El Niño','El Niño 97 (EP)','El Niño 02 (CP)']
LNnames = ['All La Niña','LN 95/96 (EP)','LN 07/08 (CP)']

fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
for i in range(len(EN)):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
    c0 = ax[i].contourf(grid.lon_rho-180,grid.lat_rho,LN[i]*mask,levels=np.linspace(-100,100,11), cmap=colmap, extend='both')
    c1 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[100],colors='k',linewidths=1)#,transform=ccrs.PlateCarree(central_longitude=180))
    c2 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[400],colors='k',linewidths=1,linestyles='dashed')#,transform=ccrs.PlateCarree(central_longitude=180))
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.3,0.1,''+str(LNnames[i])+'',transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    ax[i].set_ylim([-60,80])
    ax[i].set_extent([180-360, 295-360, -40, 40])
#    ax[i].set_xticks([190, 220, 250, 280], crs=ccrs.PlateCarree())
#    ax[i].set_yticks([-40, -20, 0, 20, 40], crs=ccrs.PlateCarree())
#    ax[i].set_ylim([-40, 40])
    ax[i].grid()
fig.subplots_adjust(right=0.85)
cb = fig.add_axes([0.9,0.25,0.01,0.5])
fig.colorbar(c0,cax=cb,label='m')
#cb.set_label('[m]',labelpad=-20, y=1.09, rotation=0)
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZisolinesatdepth_mean_LN_'+str(thresh)+'.png', bbox_inches="tight", dpi=300)

#%% Next I want to show the change in % instead of m
## Plot the mean depth of the oxyisolines for different ENSO modes
#
#EN = [rel_changeEN, rel_changeEN97, rel_changeEN02]
#LN =[rel_changeLN, rel_changeLN88, rel_changeLN07]
#ENnames = ['All El Niño','El Niño 97 (EP)','El Niño 02 (CP)']
#LNnames = ['All La Niña','La Niña 88 (EP)','La Niña 07 (CP)']
#
#fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
#for i in range(len(EN)):
#    ax[i].coastlines(color='#888888',linewidth=1.5)
#    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
#    c0 = ax[i].contourf(grid.lon_rho-180,grid.lat_rho,EN[i]*mask,levels=np.linspace(-50,50,11), cmap=cm.cm.balance, extend='both')
#    c1 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[100],colors='k',linewidths=1)
#    c2 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[400],colors='k',linewidths=1,linestyles='dashed')
#    lon_formatter = LongitudeFormatter(zero_direction_label=True)
#    lat_formatter = LatitudeFormatter()
#    ax[i].xaxis.set_major_formatter(lon_formatter)
#    ax[i].yaxis.set_major_formatter(lat_formatter)
#    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
#    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#    ax[i].text(0.7,0.92,''+str(ENnames[i])+'',transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#    ax[i].set_ylim([-60,80])
#    ax[i].set_extent([180-360, 295-360, -40, 40])
##    ax[i].set_xticks([190, 220, 250, 280], crs=ccrs.PlateCarree())
##    ax[i].set_yticks([-40, -20, 0, 20, 40], crs=ccrs.PlateCarree())
##    ax[i].set_ylim([-40, 40])
#    ax[i].grid()
#fig.subplots_adjust(right=0.85)
#cb = fig.add_axes([0.9,0.25,0.01,0.5])
#fig.colorbar(c0,cax=cb)
##cb.set_label('[m]',labelpad=-20, y=1.09, rotation=0)
##plt.subplots_adjust(wspace=0.05)
#plt.show()
##fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZisolinesatdepth_mean_EN_rel_change_'+str(thresh)+'.png', bbox_inches="tight", dpi=300)

#%% Plot to show the La Niña change of the depth of OMZ

EN = [diff_allyearsEN82, diff_allyearsEN97, diff_allyearsEN15]
ENnames = ['EN 82/83 (EP)','EN 97/98 (EP)','EN 15/16 (EP)']

fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
for i in range(len(EN)):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
    c0 = ax[i].contourf(grid.lon_rho-180,grid.lat_rho,EN[i]*mask,levels=np.linspace(-100,100,11), cmap=colmap, extend='both')
    c1 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[100],colors='k',linewidths=1)#,transform=ccrs.PlateCarree(central_longitude=180))
    c2 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[400],colors='k',linewidths=1,linestyles='dashed')#,transform=ccrs.PlateCarree(central_longitude=180)) 
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.3,0.1,''+str(ENnames[i])+'',transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    ax[i].set_ylim([-60,80])
    ax[i].set_extent([180-360, 295-360, -40, 40])
#    ax[i].set_xticks([190, 220, 250, 280], crs=ccrs.PlateCarree())
#    ax[i].set_yticks([-40, -20, 0, 20, 40], crs=ccrs.PlateCarree())
#    ax[i].set_ylim([-40, 40])
    ax[i].grid()
fig.subplots_adjust(right=0.85)
cb = fig.add_axes([0.9,0.25,0.01,0.5])
fig.colorbar(c0,cax=cb,label='m')
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZisolinesatdepth_mean_EN_EPevents'+str(thresh)+'.png', bbox_inches="tight", dpi=300)

#%% Next I want to show the three strongest EP El Niño events
# Plot the mean depth of the oxyisolines for different ENSO modes

#EN = [rel_changeEN, rel_changeEN97, rel_changeEN02]
#LN =[rel_changeLN, rel_changeLN88, rel_changeLN07]
#ENnames = ['All El Niño','El Niño 97 (EP)','El Niño 02 (CP)']
#LNnames = ['All La Niña','La Niña 88 (EP)','La Niña 07 (CP)']
#
#fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
#for i in range(len(EN)):
#    ax[i].coastlines(color='#888888',linewidth=1.5)
#    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
#    c0 = ax[i].contourf(grid.lon_rho-180,grid.lat_rho,EN[i]*mask,levels=np.linspace(-50,50,11), cmap=cm.cm.balance, extend='both')
#    c1 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[100],colors='k',linewidths=1)#,transform=ccrs.PlateCarree(central_longitude=180))
#    c2 = ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[0,:,:],levels=[400],colors='k',linewidths=1,linestyles='dashed')#,transform=ccrs.PlateCarree(central_longitude=180))
#    lon_formatter = LongitudeFormatter(zero_direction_label=True)
#    lat_formatter = LatitudeFormatter()
#    ax[i].xaxis.set_major_formatter(lon_formatter)
#    ax[i].yaxis.set_major_formatter(lat_formatter)
#    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
#    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#    ax[i].text(0.7,0.92,''+str(ENnames[i])+'',transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#    ax[i].set_ylim([-60,80])
#    ax[i].set_extent([180-360, 295-360, -40, 40])
##    ax[i].set_xticks([190, 220, 250, 280], crs=ccrs.PlateCarree())
##    ax[i].set_yticks([-40, -20, 0, 20, 40], crs=ccrs.PlateCarree())
##    ax[i].set_ylim([-40, 40])
#    ax[i].grid()
#fig.subplots_adjust(right=0.85)
#cb = fig.add_axes([0.9,0.25,0.01,0.5])
#fig.colorbar(c0,cax=cb,label='relative change (%)')
#plt.show()
##fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZisolinesatdepth_mean_EN_rel_change_'+str(thresh)+'.png', bbox_inches="tight", dpi=300)


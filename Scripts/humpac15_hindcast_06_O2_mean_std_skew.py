#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 15:04:07 2019

@author: meilers
"""

# Load packages
# Import python modules
import numpy as np
import xarray as xr
import netCDF4
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scipy.stats as sp
from scipy.interpolate import interp1d
import cmocean as cm
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
exec(open('/home/koehne/Documents/scripts/python/martinfr/myromstools.py').read())

# Define variables
depth = 100
variable = 'O2'
year_vec = range(1979,2017)

#timevector = np.zeros((1,1))
O2_timeseries = np.zeros((len(year_vec)*12,1009,701))

#
grid = xr.open_dataset('/net/kryo/work/koehne/roms/inputs/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/grd/humpac15_grd.nc')

# Open dataset
i = 0
for year in year_vec:
    data = xr.open_dataset('/net/kryo/work/koehne/roms/output/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_hindcast_1979_2016_6/monthly/his_zslice/'+str(variable)+'_'+str(depth)+'/z_humpac15_'+str(year)+'_his_'+variable+'_'+str(depth)+'.nc')
    O2_timeseries[i*12:(i+1)*12,:,:] = data.O2.values[1:,:,:]
    data.close()
    print('In year: '+str(year))
    i += 1

   
# calculate oxygen mean, std, skewness
meanO2 = O2_timeseries.mean(axis=0)
stdO2 = O2_timeseries.std(axis=0)
skewO2 = sp.stats.skew(O2_timeseries)

#%% plot mean oxygen
#fig,ax = plt.subplots()
#c0 = plt.contourf(grid.lon_rho,grid.lat_rho,meanO2,levels=np.linspace(0,300,21),extend='max')
#plt.colorbar(c0)
#plt.show()

#%% plot std oxygen
#fig,ax = plt.subplots()
#c0 = plt.contourf(grid.lon_rho,grid.lat_rho,stdO2,levels=np.linspace(0,20,11),extend='max')
#plt.colorbar(c0)
#plt.show()

#%% plot skewness oxygen
#fig,ax = plt.subplots()
#c0 = plt.contourf(grid.lon_rho,grid.lat_rho,skewO2,levels=np.linspace(-10,10,21),extend='max')
#plt.colorbar(c0)
#plt.show()

mask = np.zeros_like(grid.mask_rho) + np.NaN
mask[grid.mask_rho==1] = 1

#%% Make mask for OMZ
maskOMZ = np.zeros_like(grid.mask_rho) + np.NaN
maskOMZ[np.where(meanO2[:,:]<60)] = 1

#%%
fig,ax = plt.subplots(1,3,sharex = True, sharey = True, figsize=(12,4))
c0 = ax[0].contourf(grid.lon_rho,grid.lat_rho,meanO2,levels=np.linspace(0,300,21),extend='max')
ax[0].set_title('Depth '+str(depth)+'m, mean O2')
plt.colorbar(c0,ax=ax[0],orientation='horizontal')
c1 = ax[1].contourf(grid.lon_rho,grid.lat_rho,stdO2,levels=np.linspace(0,20,11),extend='max')
ax[1].set_title('std O2')
plt.colorbar(c1,ax=ax[1],orientation='horizontal')
c2 = ax[2].contourf(grid.lon_rho,grid.lat_rho,skewO2,levels=np.linspace(-5,5,21),extend='max',cmap=cm.cm.balance)
ax[2].set_title('skewness O2')
plt.colorbar(c2,ax=ax[2],orientation='horizontal')
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/mean_std_skew/plots/humpac15_hindcast_06_mean_std_skew/plot_mean-std-skew_'+str(depth)+'m_allpoints.png')

#%% zoom into ETP
latmin = -30
latmax = 30
lonmin=180
lonmax=300

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, sharex=True, figsize=(9,4))
#ax[0].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
#ax[1].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
#ax[2].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
#ax[0].coastlines(color='#888888',linewidth=1.5)
#ax[1].coastlines(color='#888888',linewidth=1.5)
#ax[2].coastlines(color='#888888',linewidth=1.5)
#ax[0].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
#ax[1].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
#ax[2].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))

c0 = ax[0].contourf(grid.lon_rho-180,grid.lat_rho,meanO2*mask,levels=np.linspace(0,300,21),extend='max',cmap=cm.cm.haline)
c00 = ax[0].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5,22,60],colors='k',linewidths=4)
c01 = ax[0].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5,22,60],colors='w')
plt.clabel(c00,fmt='%1d',inline_spacing=1)
plt.clabel(c01,fmt='%1d',inline_spacing=1)
ax[0].set_title('mean O2')
ax[0].set_xlim(lonmin,lonmax)
ax[0].set_ylim(latmin,latmax)
ax[0].text(265,20,'Depth: '+str(depth)+'m')
plt.colorbar(c0,ax=ax[0],orientation='horizontal')
c1 = ax[1].contourf(grid.lon_rho-180,grid.lat_rho,stdO2*mask,levels=np.linspace(0,60,11),extend='max',cmap=cm.cm.haline)
c10 = ax[1].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5,22,60],colors='k',linewidths=4)
c11 = ax[1].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5,22,60],colors='w')
plt.clabel(c10,fmt='%1d',inline_spacing=1)
plt.clabel(c11,fmt='%1d',inline_spacing=1)
ax[1].text(265,20,'Depth: '+str(depth)+'m')
ax[1].set_title('std O2')
plt.colorbar(c1,ax=ax[1],orientation='horizontal')
c2 = ax[2].contourf(grid.lon_rho-180,grid.lat_rho,skewO2*mask,levels=np.linspace(-5,5,21),extend='both',cmap=colmap)
c20 = ax[2].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5,22,60],colors='k',linewidths=4)
c21 = ax[2].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5,22,60],colors='w')
plt.clabel(c20,fmt='%1d',inline_spacing=1)
plt.clabel(c21,fmt='%1d',inline_spacing=1,colors='g')
ax[2].set_title('skewness O2')
ax[2].text(265,20,'Depth: '+str(depth)+'m')
plt.colorbar(c2,ax=ax[2],orientation='horizontal')
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/mean_std_skew/plots/humpac15_hindcast_06_mean_std_skew/plot_mean-std-skew_'+str(depth)+'m.png')

#%% New try

Var_to_plot = [meanO2,stdO2,skewO2]

fig,ax = plt.subplots(1,3, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
ax[0].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
ax[1].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
ax[2].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
for i in range(len(Var_to_plot)):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
    ax[i].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
#    ax[i].contour(grid.lon_rho-360.,grid.lat_rho,ETNP60mask,levels=[0.5],colors='w',linewidths=2,transform=ccrs.PlateCarree())
#    ax[i].contour(grid.lon_rho-360.,grid.lat_rho,ETSP60mask,levels=[0.5],colors='w',linewidths=2,transform=ccrs.PlateCarree())
#    ax[i].contour(grid.lon_rho-360.,grid.lat_rho,mask_100km,levels=[0.5],colors='w',linewidths=2,transform=ccrs.PlateCarree())
#    ax[i].text(0.6,0.65,'a)',transform=ax[i].transAxes, size=10, color='k',rotation=0.,ha="center", va="center",bbox=dict(boxstyle="square",facecolor='w',alpha=0.5))
#    ax[i].text(0.7,0.45,'b)',transform=ax[i].transAxes, size=10, color='k',rotation=0.,ha="center", va="center",bbox=dict(boxstyle="square",facecolor='w',alpha=0.5))
#    ax[i].text(0.92,0.42,'c)',transform=ax[i].transAxes, size=10, color='k',rotation=0.,ha="center", va="center",bbox=dict(boxstyle="square",facecolor='w',alpha=0.5))
    ax[i].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[60],colors='k',linewidths=1.5)
    ax[i].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[22],colors='k',linewidths=1.5,linestyles='dashed')
    ax[i].contour(grid.lon_rho-180,grid.lat_rho,meanO2*mask,[5],colors='k',linewidths=1.5,linestyles='dotted')
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([140, 180, 220, 260], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[0].text(0.7,0.92,'Mean',transform=ax[0].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    ax[1].text(0.7,0.92,'Std',transform=ax[1].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))          
    ax[2].text(0.72,0.92,'Skewness',transform=ax[2].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))        
    ax[i].set_ylim([-60,80])
    ax[i].set_extent([180-360, 295-360, -40, 40])
    ax[i].grid()
c0 = ax[0].contourf(grid.lon_rho-180,grid.lat_rho,meanO2*mask,levels=np.linspace(0,300,16), cmap=cm.cm.haline, extend='max')
c1 = ax[1].contourf(grid.lon_rho-180,grid.lat_rho,stdO2*mask,levels=np.linspace(0,60,13), cmap=cm.cm.haline, extend='max')
c2 = ax[2].contourf(grid.lon_rho-180,grid.lat_rho,skewO2*mask,levels=np.linspace(-5,5,11), cmap=colmap, extend='both')
cb0 = plt.colorbar(c0,ax=ax[0],ticks=[0,60,120,180,240,300],orientation='horizontal',label='O$_2$ (mmol m$^{-3}$)')
cb1 = plt.colorbar(c1,ax=ax[1],orientation='horizontal',label='O$_2$ (mmol m$^{-3}$)')
cb2 = plt.colorbar(c2,ax=ax[2],orientation='horizontal',label='( - )')
#ax[0].set_title('Mean O$_2$',fontsize=10)
#ax[1].set_title('Std O$_2$',fontsize=10)
#ax[2].set_title('Skewness O$_2$',fontsize=10)

plt.suptitle('a)',x=0.1,y=0.92,fontsize=12)
plt.show()
fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/mean_std_skew/plots/humpac15_hindcast_06_mean_std_skew/plot_mean-std-skew_'+str(depth)+'m_2_nocbar.png',bbox_inches="tight", dpi=300)

#%%
#fig,ax = plt.subplots(1,3, sharey=True, sharex=True, figsize=(9,4))
#c0 = ax[0].contourf(grid.lon_rho,grid.lat_rho,meanO2*mask,levels=np.linspace(0,300,21),extend='max')
#c00 = ax[0].contour(grid.lon_rho,grid.lat_rho,meanO2*mask,[5,22,60],colors='k',linewidths=4)
#c01 = ax[0].contour(grid.lon_rho,grid.lat_rho,meanO2*mask,[5,22,60],colors='w')
#plt.clabel(c00,fmt='%1d',inline_spacing=1)
#plt.clabel(c01,fmt='%1d',inline_spacing=1)
#ax[0].set_title('mean O2')
#ax[0].set_xlim(lonmin,lonmax)
#ax[0].set_ylim(latmin,latmax)
#ax[0].text(265,20,'Depth: '+str(depth)+'m')
#plt.colorbar(c0,ax=ax[0],orientation='horizontal')
#c1 = ax[1].contourf(grid.lon_rho,grid.lat_rho,stdO2*mask,levels=np.linspace(0,60,11),extend='max')
#c10 = ax[1].contour(grid.lon_rho,grid.lat_rho,meanO2*mask,[5,22,60],colors='k',linewidths=4)
#c11 = ax[1].contour(grid.lon_rho,grid.lat_rho,meanO2*mask,[5,22,60],colors='w')
#plt.clabel(c10,fmt='%1d',inline_spacing=1)
#plt.clabel(c11,fmt='%1d',inline_spacing=1)
#ax[1].text(265,20,'Depth: '+str(depth)+'m')
#ax[1].set_title('std O2')
#plt.colorbar(c1,ax=ax[1],orientation='horizontal')
#c2 = ax[2].contourf(grid.lon_rho,grid.lat_rho,skewO2*mask,levels=np.linspace(-5,5,21),extend='both',cmap=cm.cm.balance)
#c20 = ax[2].contour(grid.lon_rho,grid.lat_rho,meanO2*mask,[5,22,60],colors='k',linewidths=4)
#c21 = ax[2].contour(grid.lon_rho,grid.lat_rho,meanO2*mask,[5,22,60],colors='w')
#plt.clabel(c20,fmt='%1d',inline_spacing=1)
#plt.clabel(c21,fmt='%1d',inline_spacing=1,colors='g')
#ax[2].set_title('skewness O2')
#ax[2].text(265,20,'Depth: '+str(depth)+'m')
#plt.colorbar(c2,ax=ax[2],orientation='horizontal')
#plt.show()
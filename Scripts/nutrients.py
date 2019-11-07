# -*- coding: utf-8 -*-
"""
This is a script to visualize nutrient sources at different model depths: 10m, 50m, 100m;
includes physical, remineraliztion

Simona Meiler, October 2019, meilers@mit.edu
"""

#%% Import Python modules
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cmocean as cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

#%% Load data of the three datasets at 10m, 50m, 100m depth
data10 = xr.open_dataset('/Users/meilers/MITinternship/data/supply10m.nc')
data50 = xr.open_dataset('/Users/meilers/MITinternship/data/supply50m.nc')
data100 = xr.open_dataset('/Users/meilers/MITinternship/data/supply100m.nc')
data = xr.open_dataset('/Users/meilers/MITinternship/Data/scenario3.nc') #careful: different lon, lat!!!

lon = data10.lon    #needed for plotting
lat = data10.lat    #needed for plotting

# compile data into matrix (facilitates plotting later)
data = np.zeros((3,3,160,360))
data[0,0,:,:] = data10.transFe*31557600
data[1,0,:,:] = data50.transFe*31557600
data[2,0,:,:] = data100.transFe*31557600
data[0,1,:,:] = data10.transN*31557600
data[1,1,:,:] = data50.transN*31557600
data[2,1,:,:] = data100.transN*31557600
data[0,2,:,:] = data10.transP*31557600
data[1,2,:,:] = data50.transP*31557600
data[2,2,:,:] = data100.transP*31557600

#%% New try

# transport terms
NO3 = data50.transN
PO4 = data50.transP
FeT = data50.transFe

# other terms: remineralization and dust in case of Fe
NO3_o = data50.otherSN
PO4_o = data50.otherSP
FeT_o = data50.otherSFe

sec = 31557600  #seconds per year (365.25*86400)
#depth = 6 #RF = 0, -10, -20, -35, -55, -75, -100, -135, -185, (...)

f_dust = np.sum(data.dst_emis_solfe,axis=0)*sec/12 #in kg/m2/s --> sum up over all months, multiply with seconds of year, divide by 12 months
f_atm = f_dust/55.845 #--> convert to mol/m2/s 

N = np.add(NO3,NO3_o)*sec
P = np.add(PO4,PO4_o)*sec
Fe = np.add(FeT,FeT_o)*sec#,f_dust) #including dust

#%% Calculate ratios

# P* = Phosphat - (Nitrate/16)
P_star = data50.transP - (data50.transN/16)

# Bioavailable nutrient supply --> bio_PN, bio_FeN
# Constants
rpP = 0.0625
rpFe = 6.25e-5

bio_PN = (np.divide(P,N))*(1/rpP)
bio_FeN = (np.divide(Fe,N))*(1/rpFe)

#%% just a plot to quickly display variables

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,P,cmap=colmap)#,levels=levs[nut],cmap=colmap,extend='both')
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


#%% Plot bioavailable Fe:N and P:N

nut = 1    #chose nutrient here: 0=Fe, 1=N, 2=P
nutrient = [bio_PN,bio_FeN]
label_nut = ['N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)','Fe (mol/m$^{2}$/y)']
name_nut = ['P:N','Fe:N']

levs_PN = np.linspace(0,1000,11)
levs_FeN = np.linspace(-5,5,11)
levs = [levs_PN,levs_FeN]

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,nutrient[nut],levels=levs[nut])#,cmap=colmap,extend='both')
con1 = ax.contour(lon,lat,nutrient[nut],levels=[1.2],colors='k',linewidths=1,linstyle='solid')
con2 = ax.contour(lon,lat,nutrient[nut],levels=[2.5],colors='r',linewidths=1,linstyle='solid')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label(''+str(name_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/overview_nutr_bioav_'+str(name_nut[nut])+'.png', bbox_inches='tight', dpi=300)







#%% Plot 1 nutrient at 1 depth
nut = 2    #chose nutrient here: 0=Fe, 1=N, 2=P
depth = 1  #chose depth here: 0=10m, 1=50m, 2=100m
label_nut = ['Fe (mol/m$^{2}$/y)','N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)']
label_depth = ['10m','50m','100m']

levs_Fe = np.linspace(0,5e-09,11)
levs_N = np.linspace(0,2400,13)
levs_P = np.linspace(0,100,11)
levs = [levs_Fe,levs_N,levs_P]

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,data[depth,nut,:,:],levels=levs[nut],extend='both')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Plot 1 nutrient at 1 depth
nut = 3    #chose nutrient here: 0=Fe, 1=N, 2=P
depth = 0  #chose depth here: 0=10m, 1=50m, 2=100m
label_nut = ['Fe (mol m$^{-2}$ s$^{-1}$)','N (mol m$^{-2}$ s$^{-1}$)','P (mol m$^{-2}$ s$^{-1}$)','P* (mol m$^{-2}$ s$^{-1}$)']
label_depth = ['10m','50m','100m']

levs_Fe = np.linspace(0,5e-09,11)
levs_N = np.linspace(0,6e-05,13)
levs_P = np.linspace(0,4.5e-06,10)
levs = [levs_Fe,levs_N,levs_P]

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,P_star,levels=np.linspace(0,1e-06,11),extend='max')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label(''+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Plot nice map
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)}, sharey=True, figsize=(6,4))
#ax.imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
ax.coastlines(color='#888888',linewidth=1.5)
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
c0 = ax.contourf(lon,lat,Fe10)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.grid()
plt.colorbar(c0,ax=ax[0],orientation='vertical')
#fig.subplots_adjust(right=0.85)
#cb = fig.add_axes([0.9,0.25,0.01,0.5])
#fig.colorbar(c0,cax=cb,label='Depth (m)')
#cb.set_label('[m]',rotation=90, position=(0.5,1))
plt.show()

#%% Plot nice map - Fe across different depths
fig,ax = plt.subplots(3,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)}, sharey=True, figsize=(9,9))
#ax.imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,3):
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
c0 = ax[0].contourf(lon,lat,Fe10)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
c1 = ax[1].contourf(lon,lat,Fe50)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
c2 = ax[2].contourf(lon,lat,Fe100)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
#ax.grid()
plt.colorbar(c0,ax=ax[0],orientation='vertical')
plt.colorbar(c1,ax=ax[1],orientation='vertical')
plt.colorbar(c2,ax=ax[2],orientation='vertical')
#fig.subplots_adjust(right=0.85)

#cb = fig.add_axes([0.9,0.25,0.01,0.5])
#fig.colorbar(c0,cax=cb,label='Depth (m)')
#cb.set_label('[m]',rotation=90, position=(0.5,1))
plt.show()

#%% Plot nice map - all nutrients across same depth
fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)}, sharey=True, figsize=(9,4))
fig,ax = plt.subplots(3,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)}, sharey=True, figsize=(9,9))
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
for i in range(0,3):
    #ax[i].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    #ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
c0 = ax[0].contourf(lon,lat,Fe10)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
c1 = ax[1].contourf(lon,lat,N10)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
c2 = ax[2].contourf(lon,lat,P10)#,levels=np.logspace(-14,2,num=2), cmap=cm.cm.haline, extend='max')
#ax.grid()
plt.colorbar(c0,ax=ax[0],orientation='vertical')
plt.colorbar(c1,ax=ax[1],orientation='vertical')
plt.colorbar(c2,ax=ax[2],orientation='vertical')
plt.show()
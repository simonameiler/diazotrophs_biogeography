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

# model simulation output
Nutr = xr.open_dataset('/Users/meilers/MITinternship/Data/Nutr_tend.0000014400.nc')

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/supply50m.nc')

lon = grid.lon   #needed for plotting
lat = grid.lat    #needed for plotting

#%% 
# transport terms
NH4 = Nutr['gTr02']
NO2 = Nutr['gTr03']
NO3 = Nutr['gTr04']
PO4 = Nutr['gTr05']
FeT = Nutr['gTr07']

# other terms: remineralization and dust in case of Fe
NH4_o = Nutr['gGUD02']
NO2_o = Nutr['gGUD03']
NO3_o = Nutr['gGUD04']
PO4_o = Nutr['gGUD05']
FeT_o = Nutr['gGUD07']

#reducing the T dimension (which is 1 anyway)
NH4 = NH4[0,:,:,:]
NO2 = NO2[0,:,:,:]
NO3 = NO3[0,:,:,:]
PO4 = PO4[0,:,:,:]
FeT = FeT[0,:,:,:]

NH4_o = NH4_o[0,:,:,:]
NO2_o = NO2_o[0,:,:,:]
NO3_o = NO3_o[0,:,:,:]
PO4_o = PO4_o[0,:,:,:]
FeT_o = FeT_o[0,:,:,:]

# define constants and dz
sec = 31557600  #seconds per year (365.25*86400)
#depth = 6 #RF = 0, -10, -20, -35, -55, -75, -100, -135, -185, (...)
#dep_list = [0,10,20,35,55,75,100,135,185,260,360] #(...) and more depths
dz = [10,10,15,20,20,25] #,35,50,75,100] #(...) dz between two depth layers

#dep_m = dep_list[depth]

#%% sum nutrients up over depth and multiply with corresponding dz and sec

for i in range(len(dz)):
    NH4_int = NH4[i,:,:]*dz[i]*sec
    #print(np.min(NH4[i,:,:]))
    i+=1

for i in range(len(dz)):
    NO2_int = NO2[i,:,:]*dz[i]*sec
    i+=1
    
for i in range(len(dz)):
    NO3_int = NO3[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    PO4_int = PO4[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    FeT_int = FeT[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    NH4_int_o = NH4_o[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    NO2_int_o = NO2_o[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    NO3_int_o = NO3_o[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    PO4_int_o = PO4_o[i,:,:]*dz[i]*sec
    i+=1

for i in range(len(dz)):
    FeT_int_o = FeT_o[i,:,:]*dz[i]*sec
    i+=1

#%% Add up the different N species
N_int = NH4_int + NO3_int + NO2_int
N_int_o = NH4_int_o + NO3_int_o + NO2_int_o

#%% Set all negative values to zero

from copy import deepcopy
N_trans = deepcopy(N_int)
P_trans = deepcopy(PO4_int)
Fe_trans = deepcopy(FeT_int)
N_remin = deepcopy(N_int_o)
P_remin = deepcopy(PO4_int_o)
Fe_other = deepcopy(FeT_int_o)

#%% Set all negative values to zero - for all nutrients
# loop over depths and over latitude
for i in range(len(lat)):
    N_trans[i,:][N_int[i,:]<0]=10e-08
    #print(np.max(NH4_trans[i,:]))
    i+=1

for i in range(len(lat)):
    P_trans[i,:][PO4_int[i,:]<0]=10e-08
    i+=1

for i in range(len(lat)):
    Fe_trans[i,:][FeT_int[i,:]<0]=10e-08
    i+=1
    
for i in range(len(lat)):
    N_remin[i,:][N_int_o[i,:]<0]=10e-08
    i+=1
    
for i in range(len(lat)):
    P_remin[i,:][PO4_int_o[i,:]<0]=10e-08
    i+=1
    
for i in range(len(lat)):
    Fe_other[i,:][FeT_int_o[i,:]<0]=10e-08
    i+=1    

#%% Somehow this does not work (yet). The loop only works for the first nutrient...

#for i in range(0,10):
#    for j in range(0,22):
#        for k in range(0,159):
#            data_pos[i,j,k,:][data[i,j,k,:]<0]=0
#            print(i,j,k)
#            k+=1
#        data_pos[i,j,k,:][data[i,j,k,:]<0]=0
#        print(j)
#        j+=1
#    data_pos[i,j,k,:][data[i,j,k,:]<0]=0
#    i+=1
#    print(i)
    
#%% Define transport, remin and total nutrient fluxes; add iron flux

f_dust = np.sum(fein,axis=0)*sec # only needed to calculate the Fe remineraliztion
                                 # Fe_other contains dust. If we need Fe remin only: subtract dust
Fe_remin = Fe_other-f_dust

N_alt = np.add(N_int,N_int_o)
P_alt = np.add(PO4_int,PO4_int_o)
Fe_alt = np.add(FeT_int,FeT_int)#,f_dust) #including dust

N_tot = np.add(N_trans,N_remin)
P_tot = np.add(P_trans,P_remin)
Fe_tot = np.add(Fe_trans,Fe_other)#,f_dust) #including dust

# reminder: now we have the following variables: 
# N_tot, P_tot, Fe_tot
# N_trans, P_trans, Fe_trans 
# N_remin, P_remin, Fe_remin

#%% Calculate ratios

# Bioavailable nutrient supply --> bio_PN, bio_FeN
# Constants
rpP = 0.0625
rpFe = 6.25e-5
k = 0.1

#def div(x,y):
#    if y == 0:
#        return y
#    return x / y

bio_PN_tot = (np.divide(P_tot,N_tot))*(1/rpP)
bio_FeN_tot = (np.divide(Fe_tot,N_tot))*(1/rpFe)

bio_PN_alt = (np.divide(P_alt,N_alt))*(1/rpP)
bio_FeN_alt = (np.divide(Fe_alt,N_alt))*(1/rpFe)

bio_PN_trans = (np.divide(P_trans,N_trans))*(1/rpP)
F = np.add((Fe_trans*k),f_dust)
bio_FeN_trans = (np.divide(F,(N_trans*k)))*(1/rpFe)

bio_PN_remin = (np.divide(P_remin,N_remin))*(1/rpP)
bio_FeN_remin = (np.divide(Fe_remin,N_remin))*(1/rpFe)

#%% just a plot to quickly display variables

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,bio_FeN_alt,levels=np.linspace(0,500,11),cmap=colmap,extend='max')
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

#%% Plot 1 nutrient at 1 depth

nu = 1   #chose nutrient here: 0=Fe, 1=N, 2=P
nutr = [bio_PN_tot,bio_FeN_tot,bio_PN_trans,bio_FeN_trans,bio_PN_remin,bio_FeN_remin]
label_nut = ['N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)','Fe (mol/m$^{2}$/y)']
name_nut = ['P:N','Fe:N','transport P:N','transport Fe:N','remin P:N','remin Fe:N']

levs_PN = np.linspace(0,1000,11)
levs_FeN = np.linspace(-3,3,13)
levs = [levs_PN,levs_FeN,levs_PN,levs_FeN,levs_PN,levs_FeN]

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,np.log(nutr[nu]),levels=levs[nu],cmap=colmap,extend='both')
con1 = ax.contour(lon,lat,nutr[nu],levels=[1.2],colors='k',linewidths=1,linstyle='solid')
con2 = ax.contour(lon,lat,nutr[nu],levels=[2.5],colors='r',linewidths=1,linstyle='solid')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label(''+str(name_nut[nu])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/overview_nutr_bioav_'+str(name_nut[nut])+'.png', bbox_inches='tight', dpi=300)


#%% Plot 1 nutrient at 1 depth

nut = 2    #chose nutrient here: 0=Fe, 1=N, 2=P
nutrient = [bio_PN_tot,bio_PN_trans,bio_PN_remin]
label_nut = ['N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)','Fe (mol/m$^{2}$/y)']
name_nut = ['P:N','transport P:N','remin P:N']

levs_PN = np.linspace(0,1000,11)
levs_FeN = np.linspace(-2,3,11)
levs = [levs_PN,levs_FeN]

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,3))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
#c = ax.contourf(lon,lat,nutrient[nut])#,levels=levs[nut],cmap=colmap,extend='both')
con1 = ax.contour(lon,lat,nutr[nu],levels=[1.04],colors='purple',linewidths=1,linstyle='solid')
con2 = ax.contour(lon,lat,nutr[nut],levels=[2.5],colors='r',linewidths=1,linstyle='solid')
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
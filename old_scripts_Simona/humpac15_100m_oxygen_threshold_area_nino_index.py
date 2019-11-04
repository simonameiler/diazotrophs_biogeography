#!/usr/bin/env python

'''
Compute OMZ area for multiple thresholds and compare to ENSO

Original file by Eike Köhn, adapted by Simona Meiler, July 2019
'''
# In[1]:


# import packages
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import cartopy.crs as ccrs
import xarray as xr
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmocean as cm
import os
import sys
sys.path.append('../')
import cmocean as cm
exec(open('/home/koehne/Documents/scripts/python/martinfr/myromstools.py').read()) # import compute_zlev, compute_dz


#%% Read grid
grid_name = 'humpac15'
gridfile = '/net/kryo/work/koehne/roms/inputs/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_grd.nc'
romsGrd = getGrid(gridfile)
# Add necessary attributes
filepath = '/net/kryo/work/koehne/roms/output/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_hindcast_1979_2016_6/daily/avg/'
romsGrd.getAttrs(filepath+'humpac15_2016_avg.nc')
# Add lat lon
romsGrd.getLatLon()
# Add grid area
romsGrd.getArea()
# Add angle
romsGrd.getAngle()

grid = xr.open_dataset(gridfile)


#%%load datasets
ds = xr.open_dataset('/net/kryo/work/koehne/roms/output/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_hindcast_1979_2016_6/monthly/avg_zslice/O2_100/merged_files.nc')
print('Data Loaded')

#%%
timevec = np.arange(len(ds.time))
timevec_in_years = timevec/12.+1979

#%% Masks
mask = np.zeros_like(grid.mask_rho) + np.NaN #creates a mask to fill all the continents with NaN
mask[grid.mask_rho==1] = 1

#%%
mask_etnp = np.zeros_like(grid.lon_rho)+np.nan
lon_lims = np.logical_and(grid.lon_rho-360.<-70,grid.lon_rho-360.>-160)
lat_lims = np.logical_and(grid.lat_rho<25,grid.lat_rho>1)
mask_etnp[np.logical_and(lon_lims,lat_lims)] = 1

mask_etsp = np.zeros_like(grid.lon_rho)+np.nan
lon_lims = np.logical_and(grid.lon_rho-360.<-70,grid.lon_rho-360.>-120)
lat_lims = np.logical_and(grid.lat_rho>-25,grid.lat_rho<-1)
mask_etsp[np.logical_and(lon_lims,lat_lims)] = 1

mask_etp = np.zeros_like(grid.lon_rho)+np.nan
lon_lims = np.logical_and(grid.lon_rho-360.<-70,grid.lon_rho-360.>-160)
lat_lims = np.logical_and(grid.lat_rho>-25,grid.lat_rho<25)
mask_etp[np.logical_and(lon_lims,lat_lims)] = 1

#%% Define thresholds

thresh1 = 60
thresh2 = 22
thresh3 = 5

ex = (ds.O2<thresh1)
ex2 = (ds.O2<thresh2)
ex3 = (ds.O2<thresh3)

#%% Calculate area of OMZ

area_etnp_60 = np.zeros_like(ex)+np.nan
area_etsp_60 = np.zeros_like(ex)+np.nan
area_etp_22 = np.zeros_like(ex)+np.nan
area_etp_5 = np.zeros_like(ex)+np.nan

for i in range(np.shape(ex)[0]):
    bool_etnp_60 = ex[i,:,:]*mask_etnp*mask
    area_etnp_60[i,:,:] = bool_etnp_60*romsGrd.area
    bool_etsp_60 = ex[i,:,:]*mask_etsp*mask
    area_etsp_60[i,:,:] = bool_etsp_60*romsGrd.area    

    bool_etp_22 = ex2[i,:,:]*mask_etp*mask
    area_etp_22[i,:,:] = bool_etp_22*romsGrd.area
    
    bool_etp_5 = ex3[i,:,:]*mask_etp*mask
    area_etp_5[i,:,:] = bool_etp_5*romsGrd.area

#%%
tot_area_etnp = np.nansum(mask_etnp*romsGrd.area)
tot_area_etsp = np.nansum(mask_etsp*romsGrd.area)

tot_area_etp = np.nansum(mask_etp*romsGrd.area)

area_etnp_60_timeseries = np.zeros(len(timevec))
area_etsp_60_timeseries = np.zeros(len(timevec))
area_etp_22_timeseries = np.zeros(len(timevec))
area_etp_5_timeseries = np.zeros(len(timevec))

#%% Calculate fraction of ETP that is OMZ

for i in range(len(timevec)):
    area_etnp_60_timeseries[i] = np.nansum(area_etnp_60[i,:,:])
    area_etsp_60_timeseries[i] = np.nansum(area_etsp_60[i,:,:])
    area_etp_22_timeseries[i] = np.nansum(area_etp_22[i,:,:])
    area_etp_5_timeseries[i] = np.nansum(area_etp_5[i,:,:])


#%% Load Niño index
# read in Niño 3.4 index from model (not observed)
nino34_idx = np.load('/home/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/ENSOindices/Nino_indices/Nino34_model.npy')
# read in the normalized and smoothed (5-month running mean) Niño 3.4 index (from model)

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

#%% Plot results
    
plt.rcParams['font.size']=10
fig,ax1 = plt.subplots(figsize=(9,3.5))
#ax1.axhline(0,linestyle='-',color='#888888')
ax1.plot(timevec_in_years,(area_etsp_60_timeseries+area_etnp_60_timeseries)/(tot_area_etp),linewidth=1.5,linestyle='solid',color='C0',label='60 mmol m$^{-3}$')
ax1.plot(timevec_in_years,(area_etp_22_timeseries)/(tot_area_etp),linewidth=1,linestyle='solid',color='C1',label='22 mmol m$^{-3}$')
ax1.plot(timevec_in_years,(area_etp_5_timeseries)/(tot_area_etp),linewidth=1,linestyle='solid',color='C2',label='5 mmol m$^{-3}$')
ax1.legend(ncol=3,loc=1)
ax2 = ax1.twinx()
ax2.bar(timevec_in_years,oni_anom,width=0.1,color='#888888')
ax2.bar(timevec_in_years,oni_anom_pos,width=0.1,color='#e34a33')
ax2.bar(timevec_in_years,oni_anom_neg,width=0.1,color='#43a2ca')
ax2.axhline(0.5,linestyle='--',color='#888888')
ax2.axhline(-0.5,linestyle='--',color='#888888')
ax2.axhline(0,linestyle='-',color='#888888')
#ax1.set_xlabel('Years')
ax1.set_xlim(1979,2017)
ax2.set_ylim([-3,15])
ax1.set_ylim([-0.12,0.21])
ax2.set_yticks([-2.5,0,2.5])
ax1.set_yticks([0,0.05,0.1,0.15,0.20])  
ax2.set_ylabel('Model ONI (°C)                                           ')
ax1.set_ylabel('                      Fraction of ETP \n                     area where\n'+ r'                   O$_2$ < 60/22/5 mmol m$^{-3}$')
#ax1.grid(True,which='both')
plt.tight_layout()
plt.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZarea/plots/humpac15_hindcast_06_OMZ_area/timeseries_ETP_nino_60_22_5.png',bbox_inches='tight',dpi=300)
plt.plot()

#%% Calculate corrcoeff
testcoef60 = (area_etsp_60_timeseries+area_etnp_60_timeseries)/(tot_area_etp)
testcoef22 = (area_etp_22_timeseries)/(tot_area_etp)
testcoef5 = (area_etp_5_timeseries)/(tot_area_etp)

print(np.corrcoef(testcoef60,nino34_idx))
print(np.corrcoef(testcoef22,nino34_idx))
print(np.corrcoef(testcoef5,nino34_idx))


#%% Plot area os ETP
fig,ax = plt.subplots(1,1, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharey=True, figsize=(9,4))
    #ax[i].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
c0 = ax.contourf(grid.lon_rho-180,grid.lat_rho,mask_etp*mask,levels=np.linspace(0,1000,11), cmap=cm.cm.ice, extend='max')
#    ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[i,:,:]*mask,levels=[100],colors='w',linewidths=1)
#    ax[i].contour(grid.lon_rho-180,grid.lat_rho,data_to_plot[i,:,:]*mask,levels=[400],linestyle='dashed',colors='w',linewidths=1)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([140,160, 180,200, 220, 260], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
ax.set_yllims([-30,30])
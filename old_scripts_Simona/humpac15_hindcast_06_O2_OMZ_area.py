#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 15:14:52 2019

@author: meilers
"""

# Load packages
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
#import scipy.stats as sp
import cmocean as cm
exec(open('/home/koehne/Documents/scripts/python/martinfr/myromstools.py').read()) # import compute_zlev, compute_dz
exec(open('/home/meilers/Documents/MScThesis/Python/functions/data_load/read_nino.py').read())
exec(open('/home/meilers/Documents/MScThesis/Python/functions/data_load/read_PDO.py').read())

#%% Define variables
depths = [100,400]
thresholds = [60,22,5]
variable = 'O2'

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

#%% Read depth data (depths[0]=100m, depths[1]=400m)
depth = depths[0]
years = range(1979,2017)
O2_timeseries_100 = np.zeros((len(years)*12,1009,701))

# Open dataset and save all O2 timeseries
i = 0
for year in years:
    data = xr.open_dataset('/net/kryo/work/koehne/roms/output/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_hindcast_1979_2016_6/monthly/avg_zslice/'+str(variable)+'_'+str(depth)+'/z_humpac15_'+str(year)+'_avg_'+variable+'_'+str(depth)+'.nc')
    O2_timeseries_100[i*12:(i+1)*12,:,:] = data.O2.values[:,:,:]
    data.close()
    i += 1

#%%
depth = depths[1]
years = range(1979,2017)
O2_timeseries_400 = np.zeros((len(years)*12,1009,701))

# Open dataset and save all O2 timeseries
i = 0
for year in years:
    data = xr.open_dataset('/net/kryo/work/koehne/roms/output/humpac15/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/humpac15_hindcast_1979_2016_6/monthly/avg_zslice/'+str(variable)+'_'+str(depth)+'/z_humpac15_'+str(year)+'_avg_'+variable+'_'+str(depth)+'.nc')
    O2_timeseries_400[i*12:(i+1)*12,:,:] = data.O2.values[:,:,:]
    data.close()
    i += 1
    
#%% make mask for continents used for plotting
mask = np.zeros_like(grid.mask_rho) + np.NaN
mask[grid.mask_rho==1] = 1

#%% calculate area of OMZ
# create array for the full area and without northern part (Kamchatka) --> OMZ_area_cut
# loop over 3 thresholds and create boolean array for area below threshold
OMZ_area_full_100 = np.zeros((len(thresholds),np.shape(O2_timeseries_100)[0]))
OMZ_area_cut_100 = np.zeros_like(OMZ_area_full_100)
cut = 820
for threshidx,thresh in enumerate(thresholds):
    print(thresh)
    OMZ_bool_100 = np.where(O2_timeseries_100 < thresh, 1, 0)
    for timeidx in range(OMZ_bool_100.shape[0]):
        OMZ_area_full_100[threshidx,timeidx] = np.nansum(OMZ_bool_100[timeidx,:,:]*romsGrd.area*mask,axis=(0,1))/10**12
        OMZ_area_cut_100[threshidx,timeidx] = np.nansum(OMZ_bool_100[timeidx,:cut,:]*romsGrd.area[:cut,:]*mask[:cut,:],axis=(0,1))/10**12
        
#%%
OMZ_area_full_400 = np.zeros((len(thresholds),np.shape(O2_timeseries_400)[0]))
OMZ_area_cut_400 = np.zeros_like(OMZ_area_full_400)
cut = 820
for threshidx,thresh in enumerate(thresholds):
    print(thresh)
    OMZ_bool_400 = np.where(O2_timeseries_400 < thresh, 1, 0)
    for timeidx in range(OMZ_bool_400.shape[0]):
        OMZ_area_full_400[threshidx,timeidx] = np.nansum(OMZ_bool_400[timeidx,:,:]*romsGrd.area*mask,axis=(0,1))/10**12
        OMZ_area_cut_400[threshidx,timeidx] = np.nansum(OMZ_bool_400[timeidx,:cut,:]*romsGrd.area[:cut,:]*mask[:cut,:],axis=(0,1))/10**12
        
#%% PLOT TO SHOW WHERE THE DATA IS CUT WHEN SUMMING TO GET THE OMZ AREA
OMZ_bool_100 = np.where(O2_timeseries_400 < 60, 1, 0)
OMZ_bool_400 = np.where(O2_timeseries_400 < 60, 1, 0)

#%%
fig,ax = plt.subplots()
plt.contourf(np.mean(OMZ_bool_400,axis=0)*mask)
plt.axhline(cut, color='C1')
#plt.contourf(grid.lon_rho,grid.lat_rho,OMZ_bool[50,:,:]*mask)
#plt.hline(grid.lon_rho[cut,:],color='C1')
plt.text(0.75,0.9,'<60 mmol m$^{-3}$',transform=ax.transAxes)
ax.set_title('OMZ area '+str(depth)+'m')
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZ/plots/humpac15_hindcast_06_OMZ_area/OMZ_area_'+str(depth)+'.png')

#%%LOAD NINO3 timeseries
nino_years,nino_data = read_nino3()

#%% 
nino_idx_cho = np.where(np.logical_and(nino_years<2017,nino_years>=1979),True,False)
nino_years_cho = nino_years[nino_idx_cho]
nino_data_cho = nino_data[nino_idx_cho]

#%%LOAD PDO timeseries
PDO_years,PDO_data = read_PDO()

#%% 
PDO_idx_cho = np.where(np.logical_and(PDO_years<2017,PDO_years>=1979),True,False)
PDO_years_cho = PDO_years[PDO_idx_cho]
PDO_data_cho = PDO_data[PDO_idx_cho]

#%% PLOT THE TIMESERIES OF OMZ AREA FOR THE THREE THRESHOLDS. 

fig,ax = plt.subplots(3,1,sharex=True,figsize=(9,3))
ax[0].plot(1979+np.arange(O2_timeseries_100.shape[0])/12.,OMZ_area_full_100[0,:],color='C0',linewidth=1.5,label='60 mmol m$^{-3}$') # summation over entire pacific basin
ax[1].plot(1979+np.arange(O2_timeseries_100.shape[0])/12.,OMZ_area_cut_100[1,:],color='C1',linewidth=1.5,label='22 mmol m$^{-3}$') # summation without northern part, to discard the OMZ of Kamchatka.
ax[2].plot(1979+np.arange(O2_timeseries_100.shape[0])/12.,OMZ_area_full_100[2,:],color='C2',linewidth=1.5,label='5 mmol m$^{-3}$') # summation over entire pacific basin
for i in range(3):
    ax[i].set_xlim(1979,2017)
ax[1].set_ylabel('Area in 10$^{12}$ m$^{2}$',labelpad=10)
ax[0].set_yticks([0,4])
ax[1].set_yticks([0,1.5])
ax[2].set_yticks([0,0.2])
ax[0].set_ylim(0,8)
ax[1].set_ylim(0,3)
ax[2].set_ylim(0,0.4)
plt.suptitle('a)',x=0.08,y=0.95)
plt.subplots_adjust(hspace=0)
plt.show()
fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZarea/plots/humpac15_hindcast_06_OMZ_area/OMZ_area_100_timeline.png',bbox_inches='tight', dpi=300)

#%% PLOT THE TIMESERIES OF OMZ AREA FOR THE THREE THRESHOLDS. 

fig,ax = plt.subplots(3,1,sharex=True,figsize=(9,3))
a0, = ax[0].plot(1979+np.arange(O2_timeseries_400.shape[0])/12.,OMZ_area_cut_400[0,:],color='C0',linewidth=1.5,label='60') # summation over entire pacific basin
a1, = ax[0].plot(1979+np.arange(O2_timeseries_400.shape[0])/12.,OMZ_area_full_400[0,:],color='C0',linewidth=1.5,linestyle='dashed',label='*60') # summation over entire pacific basin
a2, = ax[1].plot(1979+np.arange(O2_timeseries_400.shape[0])/12.,OMZ_area_full_400[1,:],color='C1',linewidth=1.5,label='22') # summation without northern part, to discard the OMZ of Kamchatka.
a3, = ax[2].plot(1979+np.arange(O2_timeseries_400.shape[0])/12.,OMZ_area_full_400[2,:],color='C2',linewidth=1.5,label='5') # summation over entire pacific basin
for i in range(3):
    ax[i].set_xlim(1979,2017)
ax[1].set_ylabel('Area in 10$^{12}$ m$^{2}$',labelpad=10)
ax[0].set_yticks([36,38,40])
ax[1].set_yticks([14,16,18])
ax[2].set_yticks([3,4,5])
#ax[0].set_ylim(35,40)
#ax[1].set_ylim(10,20)
#ax[2].set_ylim(3,6)
plt.suptitle('b)',x=0.08,y=0.95)
plt.subplots_adjust(hspace=0)
ax[2].legend((a0,a1,a2,a3),('60 mmol m$^{-3}$','*60mmol m$^{-3}$','22 mmol m$^{-3}$','5 mmol m$^{-3}$'),ncol=4,loc='lower center',bbox_to_anchor=(0.5,-1))
plt.show()
fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZarea/plots/humpac15_hindcast_06_OMZ_area/OMZ_area_400_timeline.png',bbox_inches='tight', dpi=300)


#%% Plot timeseries of OMZ Area at 400m depth
fig,ax = plt.subplots(4,1,sharex=True,figsize=(10,5))
for i in range(3):
    ax[i].plot(1979+np.arange(O2_timeseries.shape[0])/12.,OMZ_area_full[i,:]) # summation over entire pacific basin
    ax[i].plot(1979+np.arange(O2_timeseries.shape[0])/12.,OMZ_area_cut[i,:]) # summation without northern part, to discard the OMZ of Kamchatka.
    ax[i].grid()
    ax[i].set_ylabel('[km$^{2}$]')
    ax[i].text(0.05,0.8,'r= '+str("{0:.2f}".format(np.corrcoef(nino_data_cho,OMZ_area_cut[i,:])[0,1])),transform=ax[i].transAxes)
    ax[i].set_title('OMZ area, '+str(depth)+'m depth, for threshold '+str(thresholds[i])+ r' mmol m$^{-3}$')
ax[-1].axhline(0,linestyle='--',color='k')
ax[-1].plot(nino_years_cho,nino_data_cho)
ax[-1].grid()
ax[-1].set_title('Nino 3 Index')
ax[-1].set_ylabel('[°C]')
ax[-1].set_xlabel('Years')
ax[-1].text(0.05,-0.7,'r = corrcoef of OMZ area and Nino3 index',transform=ax[-1].transAxes)
plt.subplots_adjust(hspace=0.4)
plt.show()
#%% PLOT THE TIMESERIES OF OMZ AREA FOR THE THREE THRESHOLDS AND WITH PDO INDEX 

fig,ax = plt.subplots(4,1,sharex=True,figsize=(10,5))
for i in range(3):
    ax[i].plot(1979+np.arange(O2_timeseries.shape[0])/12.,OMZ_area_full[i,:]) # summation over entire pacific basin
    ax[i].plot(1979+np.arange(O2_timeseries.shape[0])/12.,OMZ_area_cut[i,:]) # summation without northern part, to discard the OMZ of Kamchatka.
    ax[i].grid()
    ax[i].set_ylabel('[km$^{2}$]')
    ax[i].text(0.05,0.8,'r= '+str("{0:.2f}".format(np.corrcoef(PDO_data_cho,OMZ_area_cut[i,:])[0,1])),transform=ax[i].transAxes)
    ax[i].set_title('OMZ area, '+str(depth)+'m depth, for threshold '+str(thresholds[i])+ r' mmol m$^{-3}$')
ax[-1].axhline(0,linestyle='--',color='k')
ax[-1].plot(PDO_years_cho,PDO_data_cho)
ax[-1].grid()
ax[-1].set_title('PDO Index')
ax[-1].set_ylabel('[°C]')
ax[-1].set_xlabel('Years')
ax[-1].text(0.05,-0.7,'r = corrcoef of OMZ area and PDO index',transform=ax[-1].transAxes)
plt.subplots_adjust(hspace=0.4)
plt.show()
#fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZ/plots/humpac15_hindcast_06_OMZ_area/OMZ_area_'+str(depth)+'_timeline_PDO.png')

#%% Correlation
print(np.corrcoef(nino_data_cho,OMZ_area_cut[2,:]))
print(np.corrcoef(PDO_data_cho,OMZ_area_cut[2,:]))

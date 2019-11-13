#!/usr/bin/env python

'''
Plot OMZ Isolines 
Script to analyze spatial patterns over months --> annual variability

Compute OMZ isolines: Martin Frischknecht, Nov 2016
Plotting ideas: Eike Koehn, June 2019
Adapted by Simona Meiler, July 2019
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

#%% Run setup - needed to load saved data
setup = 'humpac15'
runtag = '{}_hindcast_1979_2016_6'.format(setup)
#inpath = '/net/kryo/work/koehne/roms/output/{}/N64ts10tb4hc250_grd_merged_SiO3_PO4_fix/{}//monthly/avg/'.format(setup,runtag)
# Set path to save output
outpath = '/net/kryo/work/meilers/Data/{}/{}/'.format(setup,runtag)
# Filelist
year_start = 1979
year_end = 2017
#filelist = [inpath+'humpac15_{}_avg.nc'.format(yr) for yr in range(year_start,year_end)]
outfile = outpath+'OMZ_depths_1979_2017_interpolated.nc'

#%%######################################################
# Load data, masks, Niño indices etc.
#########################################################

# Dataset of OMZ isolines
nc_out = netCDF4.Dataset(outfile,'r')
OMZdepth60 = nc_out['depthOMZ'][:,0,:,:]
meanOMZdepth60 = np.mean(nc_out['depthOMZ'][:,0,:,:],axis=0)

#%% read in Niño 3.4 index from model (not observed)
nino34_idx = np.load('/home/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/ENSOindices/Nino_indices/Nino34_model.npy')

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

#%% make mask for continents used for plotting
year_vec = range(1979,2017)
mask = np.zeros_like(grid.mask_rho) + np.NaN
mask[grid.mask_rho==1] = 1

area = romsGrd.area.data
depth = 100

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

fig,ax = plt.subplots()
c0 = ax.contourf(grid.lon_rho,grid.lat_rho,mask_area_cho_2D)
ax.set_xlim([225,300])
ax.set_ylim([-30,30])
ax.grid('on')

#%% Load topography for plotting
topofile='ETOPO15_mean_bath.nc'
topopath='/net/kryo/work/updata/bathymetry/ETOPO1/'
topodata = netCDF4.Dataset(topopath+topofile)
topo = topodata.variables['mean_bath'][:]
topo_lon = topodata.variables['longitude'][:]
topo_lat = topodata.variables['latitude'][:]

#%%set nc_out as omz_top to make further plotting easier

OMZdepth60 = nc_out['depthOMZ'][:,0,:,:]*-1

# make an array for means over all months
OMZ_all_months = np.zeros((12,len(grid.eta_rho),len(grid.xi_rho)))    
for i in range(12):
    OMZ_all_months[i,:,:] = np.mean(OMZdepth60[i::12,:,:],axis=0)
    
# calculate the mean over all months to later calculate anomalies from long-term mean
OMZ_mean = np.mean(OMZ_all_months, axis=0)
OMZ_mean_3D_dum = OMZ_mean[np.newaxis,:,:]
OMZ_mean_3D =np.repeat(OMZ_mean_3D_dum,12,axis=0)

#%%######################################################
# Plot
#########################################################

#levs = np.linspace(-40,40,22)
plt.rcParams['font.size'] = 10
levs = np.linspace(-300,-60,9)
fig,ax = plt.subplots(figsize=(5,4))
c0 = ax.contourf(grid.lon_rho-360.,grid.lat_rho,np.mean(OMZdepth60,axis=0),levels=levs,extend='min',cmap=cm.cm.ice_r)
        #c1 = ax[i,j].contour(grid.lon_rho-360.,grid.lat_rho,m_OMZdepth60*mask,[-350],colors='k',linestyles='-')
        #ax[i,j].set_ylim([15,30])
        #ax[i,j].set_xlim([-150,-135])
c1 = ax.contour(topo_lon,topo_lat,topo,[0],colors='k')
#c2 = ax.contour(grid.lon_rho-360.,grid.lat_rho,np.mean(OMZdepth60,axis=0),[-120,-90],colors='w')
        #ax[i,j].set_aspect('equal')
ax.set_xlim([-230,-60])
ax.set_ylim([-30,30])
ax.grid(True,which='both')
plt.colorbar(c0,ax=ax)
ax.text(-225,25, r'Mean oxycline depth in m (1979-2016)', size=10, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#plt.savefig('plots/humpac15_oxycline_average_1979-2016.png',dpi=300)


# In[11]:


#levs = np.linspace(-40,40,22)
levs = np.linspace(-300,0,11)
fig,ax = plt.subplots(3,4,figsize=(15,10))
m = 0
for i in range(3):
    for j in range(4):
        c0 = ax[i,j].contourf(grid.lon_rho-360.,grid.lat_rho,OMZ_all_months[m,:,:],levels=levs,extend='neither',cmap=cm.cm.ice_r)
        #c1 = ax[i,j].contour(grid.lon_rho-360.,grid.lat_rho,m_OMZ_all_months*mask,[-350],colors='k',linestyles='-')
        #ax[i,j].set_ylim([15,30])
        #ax[i,j].set_xlim([-150,-135])
        c1 = ax[i,j].contour(topo_lon,topo_lat,topo,[0],colors='k')
        #ax[i,j].set_aspect('equal')
        ax[i,j].set_xlim([-230,-60])
        ax[i,j].set_ylim([-30,30])
        ax[i,j].grid(True,which='both')
        m += 1    
        plt.colorbar(c0,ax=ax[i,j])


#%% Plot monthly anomalies with respect to the 1979-2016 mean
colmap = plt.get_cmap('RdBu_r')
levs = np.linspace(-40,40,21)

fig,ax = plt.subplots(3,4, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharex=True, sharey=True, figsize=(9,6))
months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
m = 0
for i in range(3):
    for j in range(4):
        #ax[i,j].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
        ax[i,j].coastlines(color='#888888',linewidth=1.5)
        ax[i,j].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
        c0 = ax[i,j].contourf(grid.lon_rho-180.,grid.lat_rho,(OMZ_all_months[m,:,:]-OMZ_mean_3D[m,:,:])*mask,levels=levs,extend='both',cmap=colmap)
        ax[i,j].contour(grid.lon_rho-180,grid.lat_rho,meanOMZdepth60[:,:],levels=[100],colors='k',linewidths=1.5)#,transform=ccrs.PlateCarree(central_longitude=180))        
        ax[i,j].contour(grid.lon_rho-180,grid.lat_rho,meanOMZdepth60[:,:],levels=[400],colors='k',linewidths=1,linestyles='dashed')#,transform=ccrs.PlateCarree(central_longitude=180))
        ax[i,j].set_extent([220-360, 295-360, -30, 30])
        ax[i,j].set_xticks([240, 270], crs=ccrs.PlateCarree())
        ax[i,j].set_yticks([-20, 0, 20], crs=ccrs.PlateCarree())
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax[2,j].xaxis.set_major_formatter(lon_formatter)
        ax[i,0].yaxis.set_major_formatter(lat_formatter)       
        ax[i,j].grid()
        ax[i,j].text(0.78,0.9, months[m], size=10, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'),transform=ax[i,j].transAxes)
        m += 1    
        #plt.suptitle('Monthly anomalies of oxycline depth (from 1979-2016 climatology)')
        fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
        cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
        cbar = fig.colorbar(c0, cax=cbar_ax)
        cbar.set_label('shallower         (m)           deeper', rotation=270,labelpad=10)
fig.savefig('/nas/meilers/Documents/MScThesis/Python/humpac15_hindcast_06_analysis/OMZisolinesatdepth/Plots/humpac15_hindcast_06/OMZdepth_monthly_anomalies_1.png',dpi=300, bbox_inches='tight')

#%% Show lon and lat line for hovmöller diagrams that follow
### Line coordinates (e.g. at latitude -5N)
#LatLine = np.ones([50])*-10. #for ETNP: -7., for ETSP: 10.
#LonLine = np.linspace(260,270,50) #for ETNP:(250,260,50) for ETSP: (268,278,50)
#LonLat_Line = np.zeros([len(LatLine),2])
#LonLat_Line[:,0] = LonLine
#LonLat_Line[:,1] = LatLine
#
## Obs Lon Lat
#lonROMS_flat = romsGrd.lon_rho.flatten()
#latROMS_flat = romsGrd.lat_rho.flatten()
#xy = np.zeros([len(latROMS_flat),2])
#xy[:,0] = lonROMS_flat
#xy[:,1] = latROMS_flat
#
#
##%% Plot monthly anomalies with respect to the 1979-2016 mean
#colmap = plt.get_cmap('RdBu_r')
#levs = np.linspace(-40,40,21)
#
#fig,ax = plt.subplots(3,4, subplot_kw={'projection':ccrs.PlateCarree(central_longitude=180)}, sharex=True, sharey=True, figsize=(9,6))
#months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
#m = 0
#for i in range(3):
#    for j in range(4):
#        ax[i,j].imshow(np.tile(np.array([[[224, 224, 224]]], dtype=np.uint8), [2, 2, 1]), origin='upper', transform=ccrs.PlateCarree(), extent=[-180, 180, -180, 180])
#        ax[i,j].coastlines(color='#888888',linewidth=1.5)
#        ax[i,j].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='face', facecolor='#808080'))
#        c0 = ax[i,j].contourf(grid.lon_rho-180.,grid.lat_rho,(OMZ_all_months[m,:,:]-OMZ_mean_3D[m,:,:])*mask,levels=levs,extend='both',cmap=colmap)
#        ax[i,j].contour(grid.lon_rho-180,grid.lat_rho,meanOMZdepth60[:,:],levels=[100],colors='k',linewidths=1)#,transform=ccrs.PlateCarree(central_longitude=180))
#        ax[i,j].plot(LonLine,LatLine,'w',linewidth=2)
#        ax[i,j].set_extent([220-360, 295-360, -30, 30])
#        ax[i,j].set_xticks([240, 270], crs=ccrs.PlateCarree())
#        ax[i,j].set_yticks([-20, 0, 20], crs=ccrs.PlateCarree())
#        lon_formatter = LongitudeFormatter(zero_direction_label=True)
#        lat_formatter = LatitudeFormatter()
#        ax[2,j].xaxis.set_major_formatter(lon_formatter)
#        ax[i,0].yaxis.set_major_formatter(lat_formatter)       
#        ax[i,j].grid()
#        ax[i,j].text(0.78,0.9, months[m], size=10, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'),transform=ax[i,j].transAxes)
#        m += 1    
#        #plt.suptitle('Monthly anomalies of oxycline depth (from 1979-2016 climatology)')
#        fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
#        cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
#        cbar = fig.colorbar(c0, cax=cbar_ax)
#        cbar.set_label('shallower         (m)           deeper', rotation=270,labelpad=10)
#
# In[14]:
#
#
#o2 = ds_dum.variables['O2'][:,:,loncho,latcho]
#t = ds_dum.variables['temp'][:,:,loncho,latcho]
#s = ds_dum.variables['salt'][:,:,loncho,latcho]
#levs_o2 = np.linspace(0,240,9)
#levs_t = np.linspace(10,30,11)
#levs_s = np.linspace(33.8,34.8,11)
#print(o2.shape)
#fig,ax = plt.subplots(1,3,figsize=(15,5),sharey=True)
#c0 = ax[0].contourf(range(12),z[:,loncho,latcho],np.transpose(o2),levels=levs_o2,extend='neither',cmap=cm.cm.haline)#-m_omz_top_3D[m,:,:]))
#c1 = ax[0].contour(range(12),z[:,loncho,latcho],np.transpose(o2),[60],colors='w')
#ax[0].clabel(c1,fmt='%2.0f')
#ax[0].set_ylim([-300,0])
#ax[0].grid(True,which='both')
#ax[0].set_xticks(np.linspace(0,11,12))
#ax[0].set_xticklabels(months,rotation=90)
#ax[0].set_ylabel('Depth in m')
#plt.colorbar(c0,ax=ax[0],orientation='horizontal',pad=0.2)
#
#c0 = ax[1].contourf(range(12),z[:,loncho,latcho],np.transpose(t),levels=levs_t,extend='neither',cmap=cm.cm.haline)#-m_omz_top_3D[m,:,:]))
#c1 = ax[1].contour(range(12),z[:,loncho,latcho],np.transpose(t),[16,28],colors='w')
#ax[1].clabel(c1,fmt='%2.0f')
#ax[1].set_ylim([-300,0])
#ax[1].grid(True,which='both')
#ax[1].set_xticks(np.linspace(0,11,12))
#ax[1].set_xticklabels(months,rotation=90)
#plt.colorbar(c0,ax=ax[1],orientation='horizontal',pad=0.2)
#
#c0 = ax[2].contourf(range(12),z[:,loncho,latcho],np.transpose(s),levels=levs_s,extend='min',cmap=cm.cm.haline)#-m_omz_top_3D[m,:,:]))
#c1 = ax[2].contour(range(12),z[:,loncho,latcho],np.transpose(s),[34.6],colors='w')
#ax[2].clabel(c1,fmt='%2.1f')
#ax[2].set_ylim([-300,0])
#ax[2].grid(True,which='both')
#ax[2].set_xticks(np.linspace(0,11,12))
#ax[2].set_xticklabels(months,rotation=90)
#plt.colorbar(c0,ax=ax[2],orientation='horizontal',pad=0.2)
#
#fig.subplots_adjust(wspace=0.15)
#ax[0].text(1,-270, 'Oxygen', size=15, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#ax[1].text(1,-270, 'Temperature', size=15, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#ax[2].text(1,-270, 'Salinity', size=15, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#
#plt.suptitle('ETNP example')
##plt.savefig('plots/humpac15_oxcline_depth_ETNP_example.png',dpi=300)


# In[15]:


#o2 = ds_dum.variables['O2'][:,:,loncho2,latcho2]
#t = ds_dum.variables['temp'][:,:,loncho2,latcho2]
#s = ds_dum.variables['salt'][:,:,loncho2,latcho2]
#levs_o2 = np.linspace(0,240,9)
#levs_t = np.linspace(10,30,11)
#levs_s = np.linspace(34.5,35.5,11)
#print(o2.shape)
#fig,ax = plt.subplots(1,3,figsize=(15,5),sharey=True)
#c0 = ax[0].contourf(range(12),z[:,loncho2,latcho2],np.transpose(o2),levels=levs_o2,extend='neither',cmap=cm.cm.haline)#-m_omz_top_3D[m,:,:]))
#c1 = ax[0].contour(range(12),z[:,loncho2,latcho2],np.transpose(o2),[60],colors='w')
#ax[0].clabel(c1,fmt='%2.0f')
#ax[0].set_ylim([-300,0])
#ax[0].grid(True,which='both')
#ax[0].set_xticks(np.linspace(0,11,12))
#ax[0].set_xticklabels(months,rotation=90)
#ax[0].set_ylabel('Depth in m')
#plt.colorbar(c0,ax=ax[0],orientation='horizontal',pad=0.2)
#
#c0 = ax[1].contourf(range(12),z[:,loncho2,latcho2],np.transpose(t),levels=levs_t,extend='neither',cmap=cm.cm.haline)#-m_omz_top_3D[m,:,:]))
#c1 = ax[1].contour(range(12),z[:,loncho2,latcho2],np.transpose(t),[16,28],colors='w')
#ax[1].clabel(c1,fmt='%2.0f')
#ax[1].set_ylim([-300,0])
#ax[1].grid(True,which='both')
#ax[1].set_xticks(np.linspace(0,11,12))
#ax[1].set_xticklabels(months,rotation=90)
#plt.colorbar(c0,ax=ax[1],orientation='horizontal',pad=0.2)
#
#c0 = ax[2].contourf(range(12),z[:,loncho2,latcho2],np.transpose(s),levels=levs_s,extend='max',cmap=cm.cm.haline)#-m_omz_top_3D[m,:,:]))
#c1 = ax[2].contour(range(12),z[:,loncho2,latcho2],np.transpose(s),[35.1],colors='w')
#ax[2].clabel(c1,fmt='%2.1f')
#ax[2].set_ylim([-300,0])
#ax[2].grid(True,which='both')
#ax[2].set_xticks(np.linspace(0,11,12))
#ax[2].set_xticklabels(months,rotation=90)
#plt.colorbar(c0,ax=ax[2],orientation='horizontal',pad=0.2)
#
#fig.subplots_adjust(wspace=0.15)
#ax[0].text(1,-270, 'Oxygen', size=15, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#ax[1].text(1,-270, 'Temperature', size=15, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#ax[2].text(1,-270, 'Salinity', size=15, rotation=0.,ha="left", va="center",bbox=dict(boxstyle="round",facecolor='w'))
#
#plt.suptitle('ETSP example')
##plt.savefig('plots/humpac15_oxcline_depth_ETSP_example.png',dpi=300)

levs = np.linspace(-.0001,.0001,22)
#levs = np.linspace(-400,0,41)
fig,ax = plt.subplots(figsize=(15,10))
m = 0
c0 = ax.contourf(grid.lon_rho-360.,grid.lat_rho,np.sum((OMZ_all_months-OMZ_mean_3D),axis=0),levels=levs,extend='neither',cmap=cm.cm.balance)
        #c1 = ax[i,j].contour(grid.lon_rho-360.,grid.lat_rho,m_omz_top*mask,[-350],colors='k',linestyles='-')
        #ax[i,j].set_ylim([15,30])
        #ax[i,j].set_xlim([-150,-135])
c1 = ax.contour(topo_lon,topo_lat,topo,[0],colors='k')
ax.set_aspect('equal')
ax.set_xlim([-230,-60])                   
plt.colorbar(c0)
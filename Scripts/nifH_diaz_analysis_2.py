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

#OJ: can read all grid info from this file
grid = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')

lon = grid.X   #needed for plotting
lat = grid.Y    #needed for plotting

area = grid.rA
dz_all = grid.drF


# Read in diazotroph data from model - same setup as for nutrients
# Simulation: run19_33
# 5 diazotroph species. according to Steph TRAC30 to TRAC34

months_vec = range(0,12)
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

    
#%% sum nutrients up over depth and multiply with corresponding dz

#OJ: shorter (xarray is smart about matching dimensions)
#diaz_int = (diaz*dz_all).sum('Z')
diaz_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz_int[:,i,:,:] = diaz[:,i,:,:]*dz_all[i].values  
    #print(np.max(diaz_int[:,i,:,:]))
    #print(i)
diaz_int = np.sum(diaz_int,axis=1)

# define std, mean and cv (to plot later)
diaz_std = deepcopy(diaz_int)
diaz_std = np.std(diaz_std,axis=0)
diaz_int = np.mean(diaz_int,axis=0)
diaz_cv = diaz_std/diaz_int

#%%############################################################################
################### Tang and Cassar database ##################################
###############################################################################

diazotroph_observations = pd.read_csv(r'/Users/meilers/MITinternship/Data/Tang_and_Cassar-2019/nifH_Gene_Integral_mod.csv')
#print(diazotroph_observations)

#OJ: why degrade precision to float32?
# Single columns of nifH database in list format
nifH_Tri = diazotroph_observations['Trichodesmium nifH Gene (x106 copies m-2)']
nifH_UCYN_A = diazotroph_observations['UCYN-A nifH Gene (x106 copies m-2)']
nifH_UCYN_B = diazotroph_observations['UCYN-B nifH Gene (x106 copies m-2)']
nifH_UCYN_C = diazotroph_observations['UCYN-C nifH Gene (x106 copies m-2)']
nifH_Richelia = diazotroph_observations['Richelia nifH Gene (x106 copies m-2)']
nifH_Calothrix = diazotroph_observations['Calothrix nifH Gene  (x106 copies m-2)']
nifH_Gamma = diazotroph_observations['Gamma nifH Gene (x106 copies/m3)']
lon_nifH = diazotroph_observations['LONGITUDE']
lat_nifH = diazotroph_observations['LATITUDE']
year = diazotroph_observations['YEAR']
month = diazotroph_observations['MONTH']

mytypes = [
    'Trichodesmium nifH Gene (x106 copies m-2)',
    'UCYN-A nifH Gene (x106 copies m-2)',
    'UCYN-B nifH Gene (x106 copies m-2)',
    'UCYN-C nifH Gene (x106 copies m-2)',
    'Richelia nifH Gene (x106 copies m-2)',
    'Calothrix nifH Gene  (x106 copies m-2)',
    'Gamma nifH Gene (x106 copies/m3)',
    ]

#OJ: this skips NaNs while summing
nifH = diazotroph_observations[mytypes]
nifH_sum = nifH.sum(1)

#%% Conversion factors for nifH abundance to biomass

tri_low = 1.2165e-05
tri_high = 3.4722e-02
UCYN_low = 1.3888e-05
UCYN_high = 3.2407e-04
Ric_low = 5.75e-07
Ric_high = 1.0916e-05

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
lon_pres = lon_nifH[presence[0]]
lat_pres = lat_nifH[presence[0]]
lon_abs = lon_nifH[absence[0]]
lat_abs = lat_nifH[absence[0]]

#%% Import lat, lon, and regions file to define ocean basins/regions
# after Teng et al., 2014

# PROBLEM 2: CREATE THE REGIONS FROM DARWIN --> SCRIPT make_mask.py
# import regions file (and lat, lon of course too) from Darwin as soon as it's ready

regions = np.fromfile('mask_darwin.int64_360x160.bin', 'int64').reshape(160, 360)
reg_lon = np.mod(lon, 360.)
reg_lat = lat

#regions = pd.read_csv('/Users/meilers/MITinternship/Data/Regions/regions.csv', header=None).values
#reg_lon = np.loadtxt('/Users/meilers/MITinternship/Data/Regions/lons.csv', delimiter=',')
#reg_lat = np.loadtxt('/Users/meilers/MITinternship/Data/Regions/lats.csv', delimiter=',')

#%% create interpolator
I = si.RegularGridInterpolator((reg_lat, reg_lon), regions, 'nearest',
                               bounds_error=False, fill_value=None)


#%% assign each observation with corresponding region

n_d = len(lon_nifH)
latlon = np.zeros((n_d, 2))
latlon[:,0] = lat_nifH
latlon[:,1] = np.mod(lon_nifH, 360.)

#repeat for different species
# 1) Trichodesmium
n_d_Tri = len(lon_nifH[nifH_Tri>0])
latlon_Tri = np.zeros((n_d_Tri, 2))
latlon_Tri[:,0] = lat_nifH[nifH_Tri>0]
latlon_Tri[:,1] = np.mod(lon_nifH[nifH_Tri>0], 360.)

# 2) UCYN_A
n_d_UCYN_A = len(lon_nifH[nifH_UCYN_A>0])
latlon_UCYN_A = np.zeros((n_d_UCYN_A, 2))
latlon_UCYN_A[:,0] = lat_nifH[nifH_UCYN_A>0]
latlon_UCYN_A[:,1] = np.mod(lon_nifH[nifH_UCYN_A>0], 360.)

# 3) UCYN_B
n_d_UCYN_B = len(lon_nifH[nifH_UCYN_B>0])
latlon_UCYN_B = np.zeros((n_d_UCYN_B, 2))
latlon_UCYN_B[:,0] = lat_nifH[nifH_UCYN_B>0]
latlon_UCYN_B[:,1] = np.mod(lon_nifH[nifH_UCYN_B>0], 360.)

# 4) Richelia
n_d_Richelia = len(lon_nifH[nifH_Richelia>0])
latlon_Richelia = np.zeros((n_d_Richelia, 2))
latlon_Richelia[:,0] = lat_nifH[nifH_Richelia>0]
latlon_Richelia[:,1] = np.mod(lon_nifH[nifH_Richelia>0], 360.)

# interpolate
nifH_reg = I(latlon).astype(int)

nifH_reg_Tri = I(latlon_Tri).astype(int)
nifH_reg_UCYN_A = I(latlon_UCYN_A).astype(int)
nifH_reg_UCYN_B = I(latlon_UCYN_B).astype(int)
nifH_reg_Richelia = I(latlon_Richelia).astype(int)

#%% stack lat, lon and assigned regions for observations into one matrix - not sure yet if that's necessary or useful...
#obs_regs = np.zeros((n_d,3))
#obs_regs[:,0] = lat_nifH
#obs_regs[:,1] = np.mod(lon_nifH, 360.)
#obs_regs[:,2] = vals_d

#%%############################################################################
################### Analyses of nifH data #####################################
###############################################################################

# PPROBLEM 3: THE ACTUAL ANALYSES
# IDEA FOR  FIRST ANALYSIS:
# show range of nifH abundances per region (boxplots) and compare it with Darwin output
# for the same region (ranges of diazotroph biomass)

# Make boxplots for the regions
regs = np.arange(0,12,1)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
for i in regs:
    ax.boxplot(nifH_sum[nifH_reg==i], labels=str(regs[i]), whis='range')


#%% What to display in the boxplots
# 1 for each region & for each species = 11 x 4 subplots
# in each subplot: 3 bars; 
#   1) range of nifH values, 2) conversion of nifH to biomass, 3) biomass from Darwin
#   for 2) show min and max of the converted nifH data. step a) pack the data of both conversion
#   factors into one list (or row of array) and display min and max (and mean) of this combined range
#   (that is to say, there will be the double of datapoints in this category.)
# formatting: make boxplots but show full range of bars (no whiskers and outliers), show mean; work
#   with two axes. left nifH, right biomass

#%% DRAFTS AND BITS OF PIECES OF CODE ... 
# try to find right index to chose regions
#np.min(obs_regs[obs_regs[:,2]==6])
#idx = obs_regs[obs_regs[:,2]==6]
#nifH_matrix[:,0][lon_nifH == obs_regs[obs_regs[:,2]==6][:,1]]

# this gives a matrix of all nifH abundances for region 6:
nifH_reg6 = nifH[nifH_reg==6]

#OJ: end
test = nifH_Tri*tri_low

lat_reg6 = obs_regs[obs_regs[:,2]==6][:,0]
lon_reg6 = obs_regs[obs_regs[:,2]==6][:,1]


#%% Presence/absence on monthly time scales
# find a way to display the data on monthly scales
pres_month = nifH_matrix[presence,3]
abs_month = nifH_matrix[absence,3]

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

#%% mask where Darwin diazotroph biomass is simulated
mask = np.where((diaz_int > 1e-04), 1, 0)
mask_out = np.where((diaz_int < 1e-04), 1, 0)

#%%############################################################################
### Quantify how many of the nifH abundances are in the predicted province ####
############### DOES NOT WORK YET !!!!!########################################

# PROBLEM 4: FIND A MEASURE FOR THE ACCURACY OF DARWIN AND CALCULATE HOW MANY
# OBSERVATIONS ARE SIMULATED CORRECTLY (IN THE PRESENCE-ABSENCE SPACE)

# careful: make sure to get lon/lat consistent!!!

#lat_corr = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
#lon_corr = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
IN = np.sum(mask[lat_pres,lon_pres])/len(lat_pres)
print(IN)

#%% Calculate accuracy for absences

OUT = np.sum(mask_out[lat_abs,lon_abs])/len(lat_abs)
print(OUT)


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

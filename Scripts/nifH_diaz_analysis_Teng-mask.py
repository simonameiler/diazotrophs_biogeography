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

#%% sum diazotrophs up over depth and multiply with corresponding dz

#OJ: shorter (xarray is smart about matching dimensions)
#diaz_int = (diaz*dz_all).sum('Z') # SM: this does not work...
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

mytypes_short = [
    'Trichodesmium nifH Gene (x106 copies m-2)',
    'UCYN-A nifH Gene (x106 copies m-2)',
    'UCYN-B nifH Gene (x106 copies m-2)',
    'Richelia nifH Gene (x106 copies m-2)',
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

conversion_low = [tri_low, UCYN_low, UCYN_low, Ric_low]
conversion_high = [tri_high, UCYN_high, UCYN_high, Ric_high]

#%% Converted nifH to biomass values
bm_Tri = nifH_Tri*tri_low

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

#regions = np.fromfile('mask_darwin.int64_360x160.bin', 'int64').reshape(160, 360)
#reg_lon = np.mod(lon, 360.)
#reg_lat = lat

regions = pd.read_csv('/Users/meilers/MITinternship/Data/Regions/regions.csv', header=None).values
reg_lon = np.loadtxt('/Users/meilers/MITinternship/Data/Regions/lons.csv', delimiter=',').astype(int)
reg_lat = np.loadtxt('/Users/meilers/MITinternship/Data/Regions/lats.csv', delimiter=',').astype(int)

# Show mask 
fig,ax = plt.subplots(figsize=(9,6))
c = ax.imshow(regions, interpolation='none')
cbar = plt.colorbar(c,ax=ax)
plt.show()

#%% create interpolator
I = si.RegularGridInterpolator((reg_lat, reg_lon), regions, 'nearest',
                               bounds_error=False, fill_value=None)


#%% assign each observation with corresponding region

n_d = len(lon_nifH)
latlon = np.zeros((n_d, 2))
latlon[:,0] = lat_nifH
latlon[:,1] = np.mod(lon_nifH, 360.)

# interpolate
nifH_reg = I(latlon).astype(int)

nifH_reg_Tri = nifH_reg[nifH_Tri>0]
nifH_reg_UCYN_A = nifH_reg[nifH_UCYN_A>0]
nifH_reg_UCYN_B = nifH_reg[nifH_UCYN_B>0]
nifH_reg_Richelia = nifH_reg[nifH_Richelia>0]


#%%############################################################################
################### Analyses of nifH data #####################################
###############################################################################

def statsfun(x, label):
    stats = {
        'med': x.median(),
        'q1': x.min(),
        'q3': x.max(),
        'whislo': x.min(),
        'whishi': x.max(),
        'mean': x.mean(),
        'label': label,
        }
    return stats

# Make boxplots for the regions
regs = np.arange(0,12,1)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
for i in regs:
    if np.sum(nifH_reg==i) > 0:
        x = nifH_sum[nifH_reg==i]
        stats = {
            'med': x.median(),
            'q1': x.min(),
            'q3': x.max(),
            'whislo': x.min(),
            'whishi': x.max(),
            'mean': x.mean(),
            'label': str(regs[i]),
            }
        ax.bxp([stats], positions=[i], showmeans=True, showfliers=False, meanline=True)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
bxpstats = []
for i in regs:
    if np.sum(nifH_reg==i) > 0:
        x = nifH_sum[nifH_reg==i]
        stats = {
            'med': x.median(),
            'q1': x.min(),
            'q3': x.max(),
            'whislo': x.min(),
            'whishi': x.max(),
            'mean': x.mean(),
            'label': str(regs[i]),
            }
        bxpstats.append(stats)

ax.bxp(bxpstats, showmeans=True, showfliers=False, meanline=True)

#%% Show results for one species in all regions
import matplotlib as mpl
mpl.rcParams['font.size'] = 10
mpl.rcParams['legend.fontsize'] = 'medium'
mpl.rcParams['figure.titlesize'] = 'medium'

# chose species (0=Trichodesmium, 1=UCYN_A, 2=UCYN_B, 3=Richlia)
species = 0
species_labels = ['Trichodesmium', 'UCYN_A', 'UCYN_B', 'Richelia']
specs_labels = ['Tri.', 'UCYN_A', 'UCYN_B', 'Richelia']
#set axes limits
#ymin = 0
#ymax = max(nifH[mytypes_short[species]])
#y2min = 0
#y2max = max(nifH[mytypes_short[species]]*tri_high)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
bxpstats = []
ax.set_ylabel('nifH Gene (x106 copies m-2)')
ax.set_title('nifH abundance and biomass: '+str(species_labels[species]))
#ax.set_ylim([ymin,ymax])
ax2 = ax.twinx()
ax2.set_ylabel('biomass (mmol C m-2)')
#ax2.set_ylim([y2min,y2max])

for i in regs:
    if np.sum(nifH_reg==i) > 0:
        nif = nifH[mytypes_short[species]][nifH_reg==i]
        bm = (nifH[mytypes_short[species]][nifH_reg==i]*tri_low).append(nifH[mytypes_short[species]][nifH_reg==i]*tri_high)
        nif_stats = statsfun(nif,str(regs[i]))
        bm_stats = statsfun(bm,str(regs[i]))
#        nif_stats = statsfun(nif,'nifH'+str(regs[i]))
#        bm_stats = statsfun(bm,'bm'+str(regs[i]))
        ax.bxp([nif_stats], positions=[i-0.25], showmeans=True, showfliers=False, meanline=True)
        ax2.bxp([bm_stats], positions=[i+0.25], showmeans=True, showfliers=False, meanline=True)

#%% Show results for all species per region

# chose regions (0-12)
reg_num = 9

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 4))
bxpstats = []
ax.set_ylabel('nifH Gene (x106 copies m-2)')
ax.set_title('nifH abundance and biomass: region '+str(reg_num))
#ax.set_ylim([ymin,ymax])
ax2 = ax.twinx()
ax2.set_ylabel('biomass (mmol C m-2)')
#ax2.set_ylim([y2min,y2max])

for i in range(0,len(mytypes_short)):
    if np.sum(nifH[mytypes_short[i]][nifH_reg==reg_num]) > 0:
        nif = nifH[mytypes_short[i]][nifH_reg==reg_num]
        bm = (nifH[mytypes_short[i]][nifH_reg==reg_num]*conversion_low[i]).append(nifH[mytypes_short[i]][nifH_reg==i]*conversion_high[i])
        nif_stats = statsfun(nif,str(specs_labels[i]))
        bm_stats = statsfun(bm,str(specs_labels[i]))
#        nif_stats = statsfun(nif,'nifH'+str(regs[i]))
#        bm_stats = statsfun(bm,'bm'+str(regs[i]))
        ax.bxp([nif_stats], positions=[i-0.25], showmeans=True, showfliers=False, meanline=True)
        ax2.bxp([bm_stats], positions=[i+0.25], showmeans=True, showfliers=False, meanline=True)

#%% show all

# FIX SECONDARY AXIS IN LOOP

fig,ax = plt.subplots(nrows=4, ncols=1,figsize=(9, 9))
bxpstats = []

for i in range(0,len(mytypes_short)):
    if np.sum(nifH[mytypes_short[i]][nifH_reg==reg_num]) > 0:
        for j in regs:
            if np.sum(nifH_reg==j) > 0:
                nif = nifH[mytypes_short[i]][nifH_reg==j]
                bm = (nifH[mytypes_short[i]][nifH_reg==j]*tri_low).append(nifH[mytypes_short[i]][nifH_reg==j]*tri_high)
                nif_stats = statsfun(nif,str(regs[j]))
                bm_stats = statsfun(bm,str(regs[j]))
#                nif_stats = statsfun(nif,'nifH'+str(regs[i]))
#                bm_stats = statsfun(bm,'bm'+str(regs[i]))
                ax[i].bxp([nif_stats], positions=[j-0.25], showmeans=True, showfliers=False, meanline=True)
                ax2 = ax[i].twinx()
                ax2.bxp([bm_stats], positions=[j+0.25], showmeans=True, showfliers=False, meanline=True)
               
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


#%%############################################################################
### Quantify how many of the nifH abundances are in the predicted province ####
###############################################################################

# mask where Darwin diazotroph biomass is simulated
mask = np.where((diaz_int > 1e-04), 1, 0)
mask_out = np.where((diaz_int < 1e-04), 1, 0)


#SM: What are the correct values for lat, lon here?
I_in = si.RegularGridInterpolator((reg_lat, reg_lon), mask, 'nearest',
                               bounds_error=False, fill_value=None)

mask_darwin_nifH = I_in(latlon).astype(int)
mask_nifH = np.where(nifH_sum > 0, 1, 0)

ncoincide = np.sum(mask_nifH == mask_darwin_nifH)
COR = ncoincide/len(nifH)
print(COR)

mask_darwin_nifH = I_in(latlon).astype(bool)
mask_nifH = nifH_sum > 0

nboth = np.sum(mask_nifH & mask_darwin_nifH)
nonlydarwin = np.sum((~mask_nifH) & mask_darwin_nifH)
nnon = np.sum((~mask_nifH) & ~mask_darwin_nifH)
nonlynifH = np.sum((mask_nifH) & ~mask_darwin_nifH)

#%% Do the math
IN = nboth/len(presence[0])
print(IN)

OUT = nnon/len(absence[0])
print(OUT)

# nifH absence but Darwin presence
False_Neg = nonlynifH/len(presence[0])
print(False_Neg)

# nifH presence but Darwin absence
False_Pos = nonlydarwin/len(absence[0])
print(False_Pos)

# check if total okay
print(IN + False_Neg)
print(OUT + False_Pos)

# SM: to calculate the values for coincide from IN and OUT - just to double check
# ((OUT*len(absence[0]))+(IN*len(presence[0])))/(len(nifH_sum))

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

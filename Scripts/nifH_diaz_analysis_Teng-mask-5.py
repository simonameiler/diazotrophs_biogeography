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

diaz1 = np.zeros((len(months_vec),23,160,360)) #smallest diazotroph --> Crocosphaera (UCYN-B)
diaz2 = np.zeros((len(months_vec),23,160,360))
diaz3 = np.zeros((len(months_vec),23,160,360))
diaz4 = np.zeros((len(months_vec),23,160,360))
diaz5 = np.zeros((len(months_vec),23,160,360)) #largest diazotroph --> Trichodesmium

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

diaz_mon = deepcopy(diaz_int)

# define std, mean and cv (to plot later)
diaz_std = deepcopy(diaz_int)
diaz_std = np.std(diaz_std,axis=0)
diaz_int = np.mean(diaz_int,axis=0)
diaz_cv = diaz_std/diaz_int

#%% Create a mean annual, depth integrated array for each diazotroph species
diaz1_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz1_int[:,i,:,:] = diaz1[:,i,:,:]*dz_all[i].values  
diaz1_int = np.sum(diaz1_int,axis=1)
diaz1_int = np.mean(diaz1_int,axis=0)

diaz2_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz2_int[:,i,:,:] = diaz2[:,i,:,:]*dz_all[i].values  
diaz2_int = np.sum(diaz2_int,axis=1)
diaz2_int = np.mean(diaz2_int,axis=0)

diaz3_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz3_int[:,i,:,:] = diaz3[:,i,:,:]*dz_all[i].values  
diaz3_int = np.sum(diaz3_int,axis=1)
diaz3_int = np.mean(diaz3_int,axis=0)

diaz4_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz4_int[:,i,:,:] = diaz4[:,i,:,:]*dz_all[i].values  
diaz4_int = np.sum(diaz4_int,axis=1)
diaz4_int = np.mean(diaz4_int,axis=0)

diaz5_int = np.zeros((12,23,160,360))
for i in range(len(dz_all)):
    diaz5_int[:,i,:,:] = diaz5[:,i,:,:]*dz_all[i].values  
diaz5_int = np.sum(diaz5_int,axis=1)
diaz5_int = np.mean(diaz5_int,axis=0)


#%% Make seasonal diazotroph biomass arrays

season = ['DJF','MAM','JJA','SON']

diaz_DJF = (diaz_mon[0,:,:]+diaz_mon[1,:,:]+diaz_mon[11,:,:])/3
diaz_MAM = np.mean(diaz_mon[2:4,:,:],axis=0)
diaz_JJA = np.mean(diaz_mon[5:7,:,:],axis=0)
diaz_SON = np.mean(diaz_mon[8:10,:,:],axis=0)

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
#nifH = diazotroph_observations[mytypes]
#nifH_sum = nifH.sum(1)

#%% Just to try and see what happens if I chose the 4 species instead of all 7 species of the dataset
# If I use this option to quantify the accuracy, I also have to adapt the var_list and exclude the 3 species there too
# Additionally, I should make sure to get the presence, absence datapoints right. They would change too.

nifH = diazotroph_observations[mytypes_short]
nifH_sum = nifH.sum(1)

# count missing values in dataframe
diazotroph_observations[mytypes_short].isnull().sum()

#%% Conversion factors for nifH abundance to biomass

tri_low = 1.2165e-05
tri_high = 3.4722e-02
UCYN_low = 1.3888e-05
UCYN_high = 3.2407e-04
Ric_low = 5.75e-07
Ric_high = 1.0916e-05

conversion_low = [tri_low, UCYN_low, UCYN_low, Ric_low]
conversion_high = [tri_high, UCYN_high, UCYN_high, Ric_high]

#%% Analyse total biomass from nifH at every point of observation

# Create empty array and fill with values of converted biomass for each species from nifH abundance
shape = len(nifH),len(mytypes_short)
bm_nifH_low = np.zeros(shape)        
bm_nifH_high = np.zeros(shape)

for i in range(0,len(mytypes_short)):
    bm_nifH_low[:,i] = nifH[mytypes_short[i]]*conversion_low[i]
bm_nifH_low_tot = np.nansum(bm_nifH_low,axis=1)

for i in range(0,len(mytypes_short)):
    bm_nifH_high[:,i] = nifH[mytypes_short[i]]*conversion_high[i]
bm_nifH_high_tot = np.nansum(bm_nifH_high,axis=1)

#%% Prepare interpolator - Regions from Darwin mask

darwin_lon = np.mod(lon, 360.)
darwin_lat = lat
#regions = np.fromfile('mask_darwin.int64_360x160.bin', 'int64').reshape(160, 360)

j, i = np.mgrid[:160, :360]

#%% create interpolator
Ii = si.RegularGridInterpolator((darwin_lat, darwin_lon), i, 'nearest',
                                bounds_error=False, fill_value=None)
Ij = si.RegularGridInterpolator((darwin_lat, darwin_lon), j, 'nearest',
                                bounds_error=False, fill_value=None)
n_d = len(lon_nifH)
latlon = np.zeros((n_d, 2))
latlon[:,0] = lat_nifH
latlon[:,1] = np.mod(lon_nifH, 360.)

# indices of darwin grid cells containing nifH observations
i_nifH = Ii(latlon).astype(int)
j_nifH = Ij(latlon).astype(int)

#%%compile the nifH data of different species into 1 matrix
#var_list = [lon_nifH, lat_nifH, year, month, nifH_Tri, nifH_UCYN_A, nifH_UCYN_B, nifH_UCYN_C, nifH_Richelia, nifH_Calothrix, nifH_Gamma]
var_list = [lon_nifH, lat_nifH, year, month, nifH_Tri, nifH_UCYN_A, nifH_UCYN_B, nifH_Richelia]
nifH_matrix = np.zeros((len(nifH_Tri),len(var_list)+1))
for i in range(len(var_list)):
    nifH_matrix[:,i] = var_list[i]

# find datapoints where no nifH gene is found at all (loop over all species)
for j in range(len(nifH_Tri)):
    if np.nansum(nifH_matrix[j,4:7]) > 0:
        nifH_matrix[j,-1] = 1

# create list of indices of nifH presences and absences. Can be used for plotting or further data manipulation
presence = np.where(nifH_matrix[:,-1] > 0)
absence = np.where(nifH_matrix[:,-1] == 0)

#%% Import lat, lon, and regions file to define ocean basins/regions
# after Teng et al., 2014

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
#latlon[:,1] = np.mod(lon_nifH, 360.)
latlon[:,1] = np.mod(lon_nifH+180.5, 360) - 180.5 # other interpolation for nifH observations as compared to the script "nifH_diaz_analysis_Teng-mask-3.py"

# interpolate
nifH_reg = I(latlon).astype(int)

nifH_reg_Tri = nifH_reg[nifH_Tri>0]
nifH_reg_UCYN_A = nifH_reg[nifH_UCYN_A>0]
nifH_reg_UCYN_B = nifH_reg[nifH_UCYN_B>0]
nifH_reg_Richelia = nifH_reg[nifH_Richelia>0]

#%% get lon, lat of presence and absence
lon_pres = latlon[:,1][presence[0]]
lat_pres = latlon[:,0][presence[0]]
lon_abs = latlon[:,1][absence[0]]
lat_abs = latlon[:,0][absence[0]]

#%% Interpolate Darwin diazotroph simulation

lm = grid.HFacC[0].values == 0

reg_lat_d = np.arange(-89., 90., 2.)
reg_lon_d = np.arange(-181., 178., 2.)

# extend for periodicity in lon
reg_lon_extended = np.r_[reg_lon_d-360, reg_lon_d, reg_lon_d+360]
regions_extended = np.c_[regions, regions, regions]

# make 2d
yr, xr = np.meshgrid(reg_lat_d, reg_lon_extended, indexing='ij')

# do not use land values
w = regions_extended != 0

# make target coordinates 2d
xo, yo = np.meshgrid(lat, (lon+182)%360-182, indexing='ij')

# map
reg_darwin = si.griddata((yr[w],xr[w]), regions_extended[w], (xo.ravel(), yo.ravel()), 'nearest')
reg_darwin = reg_darwin.reshape(160, 360)
reg_darwin[lm] = 0

#%%
fig,ax = plt.subplots(figsize=(9,6))
c = ax.imshow(reg_darwin, interpolation='none')
cbar = plt.colorbar(c,ax=ax)
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/reg_map_Darwin.png', bbox_inches='tight', dpi=100)

#%% 
#fig,ax = plt.subplots(figsize=(9,6))
#c = ax.imshow(reg_darwin==12, interpolation='none')
##cbar = plt.colorbar(c,ax=ax)
#plt.show()

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
regs = np.arange(1,13,1)

def statsfun2(x, label):
    stats = {
        'med': x.mean(),
        'q1': x.mean(),
        'q3': x.mean(),
        'whislo': x.mean(),
        'whishi': x.mean(),
        'mean': x.mean(),
        'label': label,
        }
    return stats

# chose species (0=Trichodesmium, 1=UCYN_A, 2=UCYN_B, 3=Richlia)
species = 0
species_labels = ['Trichodesmium', 'UCYN_A', 'UCYN_B', 'Richelia']
specs_labels = ['Tri.', 'UCYN_A', 'UCYN_B', 'Richelia']
#region_labels = np.arange(0,13,1)
region_labels = ['','NAtl', 'NAtlGyre', 'EqAtl', 'SAtlGyre', 'SO', 'SInd', 'NInd', 'SPacGyre', 'EqPac', 'NPacGyre', 'NPac', 'Arctic']
#set axes limits
ymin = 1e-01
ymax = 1.3e05

medianprops = dict(linestyle='-.', linewidth=0, color='k')
meanprops_Tri = dict(marker='D', markeredgecolor='black', markerfacecolor='#621055')
meanprops_A = dict(marker='D', markeredgecolor='black', markerfacecolor='#b52b65')
meanprops_B = dict(marker='D', markeredgecolor='black', markerfacecolor='#ed6663')
meanprops_Ric = dict(marker='D', markeredgecolor='black', markerfacecolor='#ffa372')
meanprops_tot = dict(marker='D', markeredgecolor='black', markerfacecolor='green')


fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 3))
bxpstats = []
ax.set_ylabel('nifH Gene (x10$^{6}$ copies m$^{-2}$)')
ax.set_title('mean nifH abundance')
ax.set_yscale('log')
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='x', bottom=False, pad=0, labelrotation=45)
#ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
ax.yaxis.grid(True, linestyle='-', which='major', color='grey',
               alpha=0.5)

pos = 0
for i in regs:
    if np.sum(nifH_reg==i) > 0:
        nif_Tri = nifH[mytypes_short[0]][nifH_reg==i]
        nif_A = nifH[mytypes_short[1]][nifH_reg==i]
        nif_B = nifH[mytypes_short[2]][nifH_reg==i]
        nif_Ric = nifH[mytypes_short[3]][nifH_reg==i]
        nif_tot = nifH_sum[nifH_reg==i]
        nif_Tri_stats = statsfun2(nif_Tri,'')
        nif_A_stats = statsfun2(nif_A,'')
        nif_B_stats = statsfun2(nif_B,str(region_labels[i]))
        nif_Ric_stats = statsfun2(nif_Ric,'')
        nif_tot_stats = statsfun2(nif_tot,'')

        c0 = ax.bxp([nif_Tri_stats], positions=[pos-0.6], showmeans=True, meanprops=meanprops_Tri, medianprops=medianprops, showfliers=False, meanline=False)
        c1 = ax.bxp([nif_A_stats], positions=[pos-0.45], showmeans=True, meanprops=meanprops_A, medianprops=medianprops, showfliers=False, meanline=False)
        c2 = ax.bxp([nif_B_stats], positions=[pos-0.3], showmeans=True, meanprops=meanprops_B, medianprops=medianprops, showfliers=False, meanline=False)
        c3 = ax.bxp([nif_Ric_stats], positions=[pos-0.15], showmeans=True, meanprops=meanprops_Ric, medianprops=medianprops, showfliers=False, meanline=False)
        c4 = ax.bxp([nif_tot_stats], positions=[pos], showmeans=True, meanprops=meanprops_tot, medianprops=medianprops, showfliers=False, meanline=False)

        ax.axvline(pos+0.2, color='k', ls='dashed',linewidth=1)
        ax.annotate(str(np.sum(nifH_reg==i)), (pos-.4,2*ymin), va='baseline', ha='center', xycoords='data')
        pos += 1

means = [c['means'][0] for c in [c0,c1,c2,c3,c4]]
ax.legend(means, 'Tri A B Ric tot'.split(),loc='center right', bbox_to_anchor=(1.12, 0.5))
#fig.savefig('/Users/meilers/MITinternship/Plots/mean_nifH_abundance_species-specific.png', bbox_inches='tight', dpi=300)
        

#%% Now plot the deducted biomass
def statsfun3(x, label):
    stats = {
        'med': x.mean(),
        'q1': x.max(),
        'q3': x.min(),
        'whislo': x.mean(),
        'whishi': x.mean(),
        'mean': x.mean(),
        'label': label,
        }
    return stats

# chose species (0=Trichodesmium, 1=UCYN_A, 2=UCYN_B, 3=Richlia)
species = 0
species_labels = ['Trichodesmium', 'UCYN_A', 'UCYN_B', 'Richelia']
specs_labels = ['Tri.', 'UCYN_A', 'UCYN_B', 'Richelia']
#set axes limits
ymin = 1e-07
ymax = 1e03

medianprops = dict(linestyle='-.', linewidth=0, color='k')

boxprops_Tri = dict(edgecolor='black', facecolor='#621055')
boxprops_A   = dict(edgecolor='black', facecolor='#b52b65')
boxprops_B   = dict(edgecolor='black', facecolor='#ed6663')
boxprops_Ric = dict(edgecolor='black', facecolor='#ffa372')
boxprops_tot = dict(edgecolor='black', facecolor='lightgreen')

bxpkw = dict(showfliers=False, showmeans=False, medianprops=medianprops, patch_artist=True)

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 3))
bxpstats = []
ax.set_ylabel('biomass (mmol C m-2)')
ax.set_title('potential biomass range from mean nifH abundance')
ax.set_yscale('log')
ax.set_ylim([ymin,ymax])
ax.tick_params(axis='x', bottom=False, pad=0, labelrotation=45)
#ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
ax.yaxis.grid(True, linestyle='-', which='major', color='grey',alpha=0.5)

#obs_count = np.zeros(len(regs))
#obs_count[0] = len(nifH_sum[nifH_reg==0])
pos = 0
for i in regs:
    if np.sum(nifH_reg==i) > 0:
        #obs_count[i] = len(nifH_sum[nifH_reg==i])
        bm_Tri = np.append(np.mean(nifH[mytypes_short[0]][nifH_reg==i])*conversion_low[0],np.mean(nifH[mytypes_short[0]][nifH_reg==i])*conversion_high[0])
        bm_A = np.append(np.mean(nifH[mytypes_short[1]][nifH_reg==i])*conversion_low[1],np.mean(nifH[mytypes_short[1]][nifH_reg==i])*conversion_high[1])
        bm_B = np.append(np.mean(nifH[mytypes_short[2]][nifH_reg==i])*conversion_low[2],np.mean(nifH[mytypes_short[2]][nifH_reg==i])*conversion_high[2])
        bm_Ric = np.append(np.mean(nifH[mytypes_short[3]][nifH_reg==i])*conversion_low[3],np.mean(nifH[mytypes_short[3]][nifH_reg==i])*conversion_high[3])
        bm_Tri_stats = statsfun3(bm_Tri,'')
        bm_A_stats = statsfun3(bm_A,'')
        bm_B_stats = statsfun3(bm_B,str(region_labels[i]))
        bm_Ric_stats = statsfun3(bm_Ric,'')

        ax.bxp([bm_Tri_stats], positions=[pos-0.6], boxprops=boxprops_Tri, **bxpkw)
        ax.bxp([bm_A_stats], positions=[pos-0.45], boxprops=boxprops_A, **bxpkw)
        ax.bxp([bm_B_stats], positions=[pos-0.3],  boxprops=boxprops_B, **bxpkw)
        ax.bxp([bm_Ric_stats], positions=[pos-0.15], boxprops=boxprops_Ric, **bxpkw)

        #ax.text(i*0.1,0.9,''+(str(obs_count[i])+''))

        ax.axvline(pos+0.2, color='k', ls='dashed',linewidth=1)
        ax.annotate(str(np.sum(nifH_reg==i)), (pos-.4,2*ymin), va='baseline', ha='center', xycoords='data')
        pos += 1

means = [c['means'][0] for c in [c0,c1,c2,c3]]
ax.legend(means, 'Tri A B Ric'.split(),loc='center right', bbox_to_anchor=(1.12, 0.5))        
#fig.savefig('/Users/meilers/MITinternship/Plots/mean_bm-from-nifH_species-specific.png', bbox_inches='tight', dpi=300)

# SM Note on results: I don't think it makes much sense to compare the mean biomass from Darwin to the biomass from nifH abundance for the 
        # different species here. Maybe first calculate a aggregate biomass estimate from nifH and then compare it to Darwin?
        # Still, these results might be biased towards the few observations...

#%% Calculate mean biomass from nifH abundance over all species 

#def statsfun4(x, label):
#    stats = {
#        'med': x.mean(),
#        'q1': x.mean(),
#        'q3': x.mean(),
#        'whislo': x.mean(),
#        'whishi': x.mean(),
#        'mean': x.mean(),
#        'label': label,
#        }
#    return stats
#
## maybe use this or a similar approach to display std around the mean of darwin biomass for each region
#def statsfun5(x, label):
#    stats = {
#        'med': x.mean(),
#        'q1': x.std(),
#        'q3': x.std(),
#        'whislo': x.mean(),
#        'whishi': x.mean(),
#        'mean': x.mean(),
#        'label': label,
#        }
#    return stats
#
#meanprops_tot_d = dict(marker='o', markeredgecolor='black', markerfacecolor='blue')
#
#ymin = 1e-03
#ymax = 1e03
#
#fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 3))
#bxpstats = []
#ax.set_ylabel('biomass (mmol C m-2)')
#ax.set_title('biomass from nifH abundance and Darwin')
#ax.set_yscale('log')
#ax.set_ylim([ymin,ymax])
#ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
#ax.yaxis.grid(True, linestyle='-', which='major', color='grey',
#               alpha=0.5)
#for i in regs:
#    if np.sum(nifH_reg==i) > 0:
#        #mean_bm_low = (np.mean(nifH[mytypes_short[0]][nifH_reg==i])*conversion_low[0])+(np.mean(nifH[mytypes_short[1]][nifH_reg==i])*conversion_low[1])+(np.mean(nifH[mytypes_short[2]][nifH_reg==i])*conversion_low[2])+(np.mean(nifH[mytypes_short[3]][nifH_reg==i])*conversion_low[3])
#        #mean_bm_high =  (np.mean(nifH[mytypes_short[0]][nifH_reg==i])*conversion_high[0])+(np.mean(nifH[mytypes_short[1]][nifH_reg==i])*conversion_high[1])+(np.mean(nifH[mytypes_short[2]][nifH_reg==i])*conversion_high[2])+(np.mean(nifH[mytypes_short[3]][nifH_reg==i])*conversion_high[3])
#        mean_bm_low = pd.concat([nifH[mytypes_short[j]][nifH_reg==i]*conversion_low[j] for j in range(4)]).mean()
#        mean_bm_high =  pd.concat([nifH[mytypes_short[j]][nifH_reg==i]*conversion_high[j] for j in range(4)]).mean()
#        mean_bm_mean = np.append(mean_bm_high,mean_bm_low)
#        bm_darwin = diaz_int[reg_darwin==i]
#        mean_bm_stats = statsfun3(mean_bm_mean,str(region_labels[i]))
#        bm_darwin_stats = statsfun4(bm_darwin,'D')
#        ax.bxp([mean_bm_stats], positions=[i-0.1], showmeans=False, patch_artist=True, medianprops=medianprops, showfliers=False, meanline=False)
#        ax.bxp([bm_darwin_stats], positions=[i+0.1], showmeans=True, meanprops=meanprops_tot_d, patch_artist=True, medianprops=medianprops, showfliers=False, meanline=False)

#fig.savefig('/Users/meilers/MITinternship/Plots/mean_bm_darwin.png', bbox_inches='tight', dpi=300)
        

#%% Analyse total biomass from nifH at every point of observation and compare to Darwin at this point

bm_darwin = diaz_int[j_nifH,i_nifH]

boxprops_darwin = dict(edgecolor='black', facecolor='lightblue')
meanprops_bm = dict(marker='D', markeredgecolor='black', markerfacecolor='green')
meanprops_dar = dict(marker='D', markeredgecolor='black', markerfacecolor='blue')
bxpkw2 = dict(showfliers=False, showmeans=True, meanprops=meanprops_bm, medianprops=medianprops, patch_artist=True)
bxpkw2dar = dict(showfliers=False, showmeans=True, meanprops=meanprops_dar, medianprops=medianprops, patch_artist=True)

ymin = 1e-04
ymax = 1e04

fig,ax = plt.subplots(nrows=1, ncols=1, figsize=(9, 3))
bxpstats = []
ax.set_ylabel('biomass (mmol C m-2)')
ax.set_title('biomass from nifH abundance and Darwin')
ax.set_yscale('log')
ax.set_ylim([ymin,ymax])
#ax.xaxis.grid(True, linestyle='-', which='major', color='lightgrey',
#               alpha=0.5)
ax.yaxis.grid(True, linestyle='-', which='major', color='grey',
               alpha=0.5)
ax.tick_params(axis='x', bottom=False, pad=0, labelrotation=45)

pos = 0
for i in regs:
    if np.sum(nifH_reg==i) > 0:
        mean_bm_low = bm_nifH_low_tot[nifH_reg==i]
        mean_bm_high =  bm_nifH_high_tot[nifH_reg==i]
        mean_bm_mean = np.append(mean_bm_high,mean_bm_low)
        mean_bm_darwin = bm_darwin[nifH_reg==i]
        mean_bm_stats = statsfun3(mean_bm_mean,str(region_labels[i]))
        bm_darwin_stats = statsfun3(mean_bm_darwin,'')
        bm0 = ax.bxp([mean_bm_stats], positions=[pos-0.1], boxprops=boxprops_tot, **bxpkw2)
        bm1 = ax.bxp([bm_darwin_stats], positions=[pos+0.1], boxprops=boxprops_darwin, **bxpkw2dar)

        ax.axvline(pos+0.5, color='k', ls='dashed',linewidth=1)
        ax.annotate(str(np.sum(nifH_reg==i)), (pos+.2,2*1e03), va='baseline', ha='center', xycoords='data')
        pos += 1

means = [c['means'][0] for c in [bm0,bm1]]
ax.legend(means, 'nifH model'.split(),loc='center right', bbox_to_anchor=(1.15, 0.5))

#fig.savefig('/Users/meilers/MITinternship/Plots/bm_darwin_nifH.png', bbox_inches='tight', dpi=300)


#%% Presence/absence on monthly time scales
# find a way to display the data on monthly scales
pres_month = nifH_matrix[presence,3]
abs_month = nifH_matrix[absence,3]

#%% Manipulating the nifH data to bring it into mappable form - SPECIES SPECIFIC

#Tri_list = np.where(nifH_Tri > 0)
#UCYN_A_list = np.where(nifH_UCYN_A > 0) 
#UCYN_B_list = np.where(nifH_UCYN_B > 0)
#UCYN_C_list = np.where(nifH_UCYN_C > 0)
#Richelia_list = np.where(nifH_Richelia > 0)
#Calothrix_list = np.where(nifH_Calothrix > 0)
#Gamma_list = np.where(nifH_Gamma > 0)
#
#
##%% Absence of nifH (zeros or n.d.)
#
#no_Tri_list = np.where(nifH_Tri == 0)
#no_UCYN_A_list = np.where(nifH_UCYN_A == 0) 
#no_UCYN_B_list = np.where(nifH_UCYN_B == 0)
#no_UCYN_C_list = np.where(nifH_UCYN_C == 0)
#no_Richelia_list = np.where(nifH_Richelia == 0)
#no_Calothrix_list = np.where(nifH_Calothrix == 0)
#no_Gamma_list = np.where(nifH_Gamma == 0)


#%%############################################################################
### Quantify how many of the nifH abundances are in the predicted province ####
###############################################################################

# mask where Darwin diazotroph biomass is simulated
thresh = 1e-05

mask = np.where((diaz_int > thresh), 1, 0)
mask_out = np.where((diaz_int < thresh), 1, 0)


#SM: What are the correct values for lat, lon here?
I_in = si.RegularGridInterpolator((lat, lon), mask, 'nearest',
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


#%% Repeat the presence/absence comparison on seasonal timescales

# reminder: diaz_DJF, diaz_JJA, diaz_MAM, diaz_SON are the seasonal arrays for the model output
seas_handle = 3 # chose season here

diaz_seasonal = [diaz_DJF, diaz_JJA, diaz_MAM, diaz_SON]

# note, months in the dataset are reported from 1-12 (not 0-11)
DJF = (diazotroph_observations['MONTH']==12) | (diazotroph_observations['MONTH']==1) | (diazotroph_observations['MONTH']==2)
JJA = (diazotroph_observations['MONTH']==3) | (diazotroph_observations['MONTH']==4) | (diazotroph_observations['MONTH']==5)
MAM = (diazotroph_observations['MONTH']==6) | (diazotroph_observations['MONTH']==7) | (diazotroph_observations['MONTH']==8)
SON = (diazotroph_observations['MONTH']==9) | (diazotroph_observations['MONTH']==10) | (diazotroph_observations['MONTH']==11)

seasons_list = [DJF, JJA, MAM, SON]

nifH_seas = diazotroph_observations[mytypes_short][seasons_list[seas_handle]]
nifH_sum_seas = nifH_seas.sum(1)

mask_seas = np.where((diaz_seasonal[seas_handle] > thresh), 1, 0)

I_in_seas = si.RegularGridInterpolator((lat, lon), mask_seas, 'nearest',
                               bounds_error=False, fill_value=None)

mask_darwin_nifH_seas = I_in_seas(latlon).astype(int)
mask_nifH_seas = np.where(seasons_list[seas_handle] & nifH_sum > 0, 1, 0)

ncoincide_seas = np.sum(mask_nifH_seas == mask_darwin_nifH_seas)
COR = ncoincide_seas/len(nifH)
print(COR)

# presence/absence on seasonal scales to plot
presence_seas = np.where(mask_nifH_seas)
abs_nifH = np.where(nifH_matrix[:,-1] == 0,1,0)
mask_nifH_seas_abs = np.where(seasons_list[seas_handle] & abs_nifH, 1, 0)
absence_seas = np.where(mask_nifH_seas_abs)

#%% Repeat the presence/absence comparison on seasonal timescales

# reminder: diaz_DJF, diaz_JJA, diaz_MAM, diaz_SON are the seasonal arrays for the model output
diaz_seasonal = [diaz_DJF, diaz_JJA, diaz_MAM, diaz_SON]

# note, months in the dataset are reported from 1-12 (not 0-11)
DJF = (diazotroph_observations['MONTH']==12) | (diazotroph_observations['MONTH']==1) | (diazotroph_observations['MONTH']==2)
JJA = (diazotroph_observations['MONTH']==3) | (diazotroph_observations['MONTH']==4) | (diazotroph_observations['MONTH']==5)
MAM = (diazotroph_observations['MONTH']==6) | (diazotroph_observations['MONTH']==7) | (diazotroph_observations['MONTH']==8)
SON = (diazotroph_observations['MONTH']==9) | (diazotroph_observations['MONTH']==10) | (diazotroph_observations['MONTH']==11)

seasons_list = [DJF, JJA, MAM, SON]

plt.rcParams.update({'font.size': 10})
col = cm.cm.haline
pres_cols = ['red','orange','purple','magenta']
abs_cols = ['lightgreen','green','lime','w']


fig,ax = plt.subplots(4,1,subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(9,12),sharex=True,sharey=True)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()

for i in range(0,4):
    nifH_seas = diazotroph_observations[mytypes_short][seasons_list[i]]
    nifH_sum_seas = nifH_seas.sum(1)

    mask_seas = np.where((diaz_seasonal[i] > thresh), 1, 0)
    mask_darwin_nifH_seas = I_in_seas(latlon).astype(int)
    mask_nifH_seas = np.where(seasons_list[i] & nifH_sum > 0, 1, 0)

    ncoincide_seas = np.sum(mask_nifH_seas == mask_darwin_nifH_seas)
    COR = ncoincide_seas/len(nifH)
    #print(COR)

# presence/absence on seasonal scales to plot
    presence_seas = np.where(mask_nifH_seas)
    abs_nifH = np.where(nifH_matrix[:,-1] == 0,1,0)
    mask_nifH_seas_abs = np.where(seasons_list[i] & abs_nifH, 1, 0)
    absence_seas = np.where(mask_nifH_seas_abs)

    ax[i].coastlines(color='#888888',linewidth=1.5)
    ax[i].add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
    ax[i].xaxis.set_major_formatter(lon_formatter)
    ax[i].yaxis.set_major_formatter(lat_formatter)
    ax[i].set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
    ax[i].set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
    ax[i].text(0.8,0.9,''+(str(season[i])+''),transform=ax[i].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
    ax[i].text(1.05,0.5,'coincidence: '+(str("{:.2%}".format(COR)+'')),transform=ax[i].transAxes, size=10, rotation=90.,ha="center", va="center")#,bbox=dict(boxstyle="square",facecolor='w'))
    ax[i].contourf(lon,lat,mask_seas,cmap=col,extend='max')
    ax[i].plot(lon_nifH[presence_seas[0]],lat_nifH[presence_seas[0]],'.',color='orange',label='nifH present')#label='Trichodesmium')
    ax[i].plot(lon_nifH[absence_seas[0]],lat_nifH[absence_seas[0]],'x',color='m',label='nifH below LOD')
    ax[i].legend(loc='lower center',ncol=2)#,bbox_to_anchor=(1.15, 1.0))

#ax[0].set_title('Seasonal nutrient ratio P:N')
#cbar.set_label(''+str(name_nut[nu])+'',rotation=90, position=(0.5,0.5))
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_mask-in-out_seas.png', bbox_inches='tight', dpi=300)

#%% calculate the percentage of coincidence as function of the threshold

#Fe:N
new_thresh = np.linspace(1e-5,1e-1,5) # choose the new ratio for the comparison here
new_thresh = [1e-5,1e-4,1e-3,1e-2,1e-1]
new_corr = np.zeros_like(new_thresh)

for i in range(len(new_thresh)):
    mask = np.where((diaz_int > new_thresh[i]), 1, 0)
    I_in = si.RegularGridInterpolator((lat, lon), mask, 'nearest', bounds_error=False, fill_value=None)
    mask_darwin_nifH = I_in(latlon).astype(int)
    mask_nifH = np.where(nifH_sum > 0, 1, 0)

    ncoincide = np.sum(mask_nifH == mask_darwin_nifH)
    new_corr[i] = ncoincide/len(nifH)


#%% plot the percentage of coincidence as function of the threshold

fig,ax = plt.subplots(1,1,figsize=(4,4),sharey=True)
ax.plot(new_corr, new_thresh)#,levels=np.linspace(0.5e14,3.5e14,7),extend='both')
#ax.axhline(1e-4,linewidth=1.0,linestyle='dashed',color='k')
#ax.axvline(ref_PN,linewidth=1.0,linestyle='dashed',color='w')
ax.set_xlabel('coincidence')
ax.set_ylabel('threshold')
ax.set_yscale('log')
#ax.text(0.85,0.95,'accuracy',transform=ax[1].transAxes, size=10, rotation=0.,ha="center", va="center",bbox=dict(boxstyle="round",facecolor='w'))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/thresh-vs-coin.png', bbox_inches='tight', dpi=300)

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
ax.plot(lon_nifH[presence[0]],lat_nifH[presence[0]],'.',color='orange',label='nifH present')#label='Trichodesmium')
ax.plot(lon_nifH[absence[0]],lat_nifH[absence[0]],'x',color='m',label='nifH below LOD')
ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifH_4species.png', bbox_inches='tight', dpi=300)

#%%

#col = plt.get_cmap('RdBu_r')
col = cm.cm.haline
from matplotlib.colors import LogNorm

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
c0 = ax.contourf(lon,lat,reg_darwin,levels=np.linspace(0,12,13),cmap=cm.cm.gray,extend='max',alpha=0.7)
c1 = ax.scatter(lon_nifH,lat_nifH,s=10,c=nifH_sum,norm=plt.Normalize(0,10000),cmap=col)
#c1 = ax.scatter(lon_nifH,lat_nifH,s=10,c=nifH_sum,norm=LogNorm(),cmap=col)
#ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c1, cax=cbar_ax)
cbar.set_label('nifH gene abundance (x10$^{6}$ copies m$^{-2}$)',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifH_regions.png', bbox_inches='tight', dpi=300)


#%% 

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
c0 = ax.contourf(lon,lat,mask,cmap=col,extend='max')
ax.plot(lon_nifH[presence[0]],lat_nifH[presence[0]],'.',color='orange',label='nifH present')#label='Trichodesmium')
ax.plot(lon_nifH[absence[0]],lat_nifH[absence[0]],'x',color='m',label='nifH below LOD')
ax.legend(loc='best')
#fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
#cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
#cbar = fig.colorbar(c0, cax=cbar_ax)
#cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_mask-in-out.png', bbox_inches='tight', dpi=300)


#%% Plot the seasonal presence/absence analysis

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
c0 = ax.contourf(lon,lat,mask_seas,cmap=col,extend='max')
ax.plot(lon_nifH[presence_seas[0]],lat_nifH[presence_seas[0]],'.',color='orange',label='nifH present')#label='Trichodesmium')
ax.plot(lon_nifH[absence_seas[0]],lat_nifH[absence_seas[0]],'x',color='m',label='nifH below LOD')
ax.legend(loc='best')
#fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
#cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
#cbar = fig.colorbar(c0, cax=cbar_ax)
#cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_mask-in-out.png', bbox_inches='tight', dpi=300)

#%%

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
c0 = ax.contourf(lon,lat,diaz_int,levels=np.linspace(0,40,21),cmap=cm.cm.gray_r,extend='max',alpha=1)
#ax.plot(lon_nifH[nifH_sum],lat_nifH[nifH_sum],'.',color='orange',label='any nifH present')#label='Trichodesmium')
ax.scatter(lon_nifH,lat_nifH,s=10,c=nifH_sum,norm=plt.Normalize(0,40),cmap=col)#label='Trichodesmium')
#ax.plot(lon_nifH[absence[0]],lat_nifH[absence[0]],'x',color='m',label='all nifH absent')
#ax.legend(loc='best')
fig.subplots_adjust(wspace=0.07,hspace=0.07,right=0.85)
cbar_ax = fig.add_axes([0.87, 0.12, 0.02, 0.75])
cbar = fig.colorbar(c0, cax=cbar_ax)
cbar.set_label('mmolC m$^{-2}$',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/diaz_Darwin_overview_nifH_regions.png', bbox_inches='tight', dpi=300)

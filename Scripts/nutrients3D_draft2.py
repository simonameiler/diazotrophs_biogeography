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
area_info = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')

lon = grid.lon   #needed for plotting
lat = grid.lat    #needed for plotting

area = area_info.rA

#%% Load and condense diazotroph data
ds = xr.open_dataset('/Users/meilers/MITinternship/Data/MarEDat20130403Diazotrophs.nc',decode_times=False)

# extract variables which are needed and convert/integrate
lon_d = ds['LONGITUDE']
lat_d = ds['LATITUDE'][10:-10] # to match the latitude of the nutrient data

obs = ds['OBSERVATIONS']
abund = ds['ABUNDANCE']
bm = ds['BIOMASS']
nifH = ds['nifHbiom']
nz_obs = ds['NON_ZERO_OBS']
nz_abund = ds['NON_ZERO_ABUND']
nz_bm = ds['NON_ZERO_BIOM']
nz_nifH = ds['NON_ZERO_nifH']

obs_tot = np.sum(obs[:,0:6,10:-10,:],axis=(0,1))
abund_tot = np.sum(abund[:,0:6,10:-10,:],axis=(0,1))
bm_tot = np.sum(bm[:,0:6,10:-10,:],axis=(0,1))
nifH_tot = np.sum(nifH[:,0:6,10:-10,:],axis=(0,1))
nz_obs_tot = np.sum(nz_obs[:,0:6,10:-10,:],axis=(0,1))
nz_abund_tot = np.sum(nz_abund[:,0:6,10:-10,:],axis=(0,1))
nz_bm_tot = np.sum(nz_bm[:,0:6,10:-10,:],axis=(0,1))
nz_nifH_tot = np.sum(nz_nifH[:,0:6,10:-10,:],axis=(0,1))

#%% Load and condense nutrient data
# transport terms
NH4 = Nutr['gTr02']
NO2 = Nutr['gTr03']
NO3 = Nutr['gTr04']
PO4 = Nutr['gTr05']
FeT = Nutr['gTr07']

# other terms: remineralization and dust in case of Fe
S_DIN = Nutr['S_DIN']
S_PO4 = Nutr['S_PO4']
S_Fe = Nutr['S_Fe']

#reducing the T dimension (which is 1 anyway)
NH4 = NH4[0,:,:,:]
NO2 = NO2[0,:,:,:]
NO3 = NO3[0,:,:,:]
PO4 = PO4[0,:,:,:]
FeT = FeT[0,:,:,:]

S_DIN = S_DIN[0,:,:,:]
S_PO4 = S_PO4[0,:,:,:]
S_Fe = S_Fe[0,:,:,:]

# define constants and dz
sec = 31557600  #seconds per year (365.25*86400)
dz = [10,10,15,20,20,25] #,35,50,75,100] #(...) dz between two depth layers

#%% sum nutrients up over depth and multiply with corresponding dz and sec

NH4_int = np.zeros((6,160,360))
for i in range(len(dz)):
    NH4_int[i,:,:] = NH4[i,:,:]*dz[i]*sec
#    print(np.max(NH4_int[i,:,:]))
#    print(i)
NH4_int = np.sum(NH4_int,axis=0)

NO2_int = np.zeros((6,160,360))
for i in range(len(dz)):
    NO2_int[i,:,:] = NO2[i,:,:]*dz[i]*sec    
NO2_int = np.sum(NO2_int,axis=0)
 
NO3_int = np.zeros((6,160,360))
for i in range(len(dz)):
    NO3_int[i,:,:] = NO3[i,:,:]*dz[i]*sec    
NO3_int = np.sum(NO3_int,axis=0)

PO4_int = np.zeros((6,160,360))
for i in range(len(dz)):
    PO4_int[i,:,:] = PO4[i,:,:]*dz[i]*sec    
PO4_int = np.sum(PO4_int,axis=0)

FeT_int = np.zeros((6,160,360))
for i in range(len(dz)):
    FeT_int[i,:,:] = FeT[i,:,:]*dz[i]*sec    
FeT_int = np.sum(FeT_int,axis=0)

S_DIN_int = np.zeros((6,160,360))
for i in range(len(dz)):
    S_DIN_int[i,:,:] = S_DIN[i,:,:]*dz[i]*sec    
S_DIN_int = np.sum(S_DIN_int,axis=0)

S_PO4_int = np.zeros((6,160,360))
for i in range(len(dz)):
    S_PO4_int[i,:,:] = S_PO4[i,:,:]*dz[i]*sec    
S_PO4_int = np.sum(S_PO4_int,axis=0)

S_Fe_int = np.zeros((6,160,360))
for i in range(len(dz)):
    S_Fe_int[i,:,:] = S_Fe[i,:,:]*dz[i]*sec    
S_Fe_int = np.sum(S_Fe_int,axis=0)

#%% Add up the different N species
N_int = NH4_int + NO3_int + NO2_int
N_int_o = S_DIN_int

#%% Set all negative values to zero

from copy import deepcopy
N_trans = deepcopy(N_int)
P_trans = deepcopy(PO4_int)
Fe_trans = deepcopy(FeT_int)
N_remin = deepcopy(N_int_o)
P_remin = deepcopy(S_PO4_int)
Fe_other = deepcopy(S_Fe_int)
  
#%% Define transport, remin and total nutrient fluxes; add iron flux

f_dust = np.sum(fein,axis=0)*sec # only needed to calculate the Fe remineraliztion
                                 # Fe_other contains dust. If we need Fe remin only: subtract dust
Fe_remin = Fe_other-f_dust

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
ref_PN = 1.04
ref_FeN = 1.2

bio_PN_tot = (np.divide(P_tot,N_tot))*(1/rpP)
bio_FeN_tot = (np.divide(Fe_tot,N_tot))*(1/rpFe)

bio_PN_trans = (np.divide(P_trans,N_trans))*(1/rpP)
F = np.add((Fe_trans*k),f_dust*100)
bio_FeN_trans = (np.divide(F,(N_trans*k)))*(1/rpFe)

bio_PN_remin = (np.divide(P_remin,N_remin))*(1/rpP)
bio_FeN_remin = (np.divide(Fe_remin,N_remin))*(1/rpFe)

#%% mask where P:N OR Fe:N is sufficient to support diazotrophs
mask = np.where((bio_FeN_tot[:,:] > ref_FeN) & (bio_PN_tot[:,:] > ref_PN), 1, 0)

#%% Manipulate diazotroph data

# Create a mask of the provinces where diazotrophs are predicted from nutrient ratios        
diaz_obs = np.zeros_like(obs_tot)
diaz_obs[obs_tot>0] = 1

diaz_abund = np.zeros_like(abund_tot)
diaz_abund[abund_tot>0] = 1

diaz_bm = np.zeros_like(bm_tot)
diaz_bm[bm_tot>0] = 1

diaz_nifH = np.zeros_like(nifH_tot)
diaz_nifH[nifH_tot>0] = 1

diaz_nz_obs = np.zeros_like(nz_obs_tot)
diaz_nz_obs[nz_obs_tot>0] = 1

diaz_nz_abund = np.zeros_like(nz_abund_tot)
diaz_nz_abund[nz_abund_tot>0] = 1

diaz_nz_bm = np.zeros_like(nz_bm_tot)
diaz_nz_bm[nz_bm_tot>0] = 1

diaz_nz_nifH = np.zeros_like(nz_nifH_tot)
diaz_nz_nifH[nz_nifH_tot>0] = 1

find_obs = np.where(diaz_obs==1)
find_abund = np.where(diaz_abund==1)
find_bm = np.where(diaz_bm==1)
find_nifH = np.where(diaz_nifH==1)
find_nz_obs = np.where(diaz_nz_obs==1)
find_nz_abund = np.where(diaz_nz_abund==1)
find_nz_bm = np.where(diaz_nz_bm==1)
find_nz_nifH = np.where(diaz_nz_nifH==1)

#pack the masks for the different diazotroph variables into one list
diaz_data_list = [find_obs,find_abund,find_bm,find_nifH,find_nz_obs,find_nz_abund,find_nz_bm,find_nz_nifH]

#%% Quantify how many of the diazotrophs abundances are in the predicted province
# careful: make sure to get lon/lat of nutrients and diazotrophs consistent!!!
# correct the two scales of latitude to match one another. (lon would be the same but to avoid confusion
# I converted it too.) All we care about here is getting the right indices matching the lon, lat of both,
# diazotroph and nutrient data. 
list_idx = 1 #to chose which data from diaz_data_list to plot

lat_corr = diaz_data_list[list_idx][0]
#lon_corr = diaz_data_list[list_idx][1]-180 #would also work...the option with %360 is nicer though
lon_corr = (diaz_data_list[list_idx][1]-180)%360
# gives fraction of abundances that are within the predicted province
IN = np.sum(mask[lat_corr,lon_corr])/len(lat_corr)
print(IN)

#%% Calculate differences in area for P:N of 0.99 to 1.04 - or any other value
new_PN = 1.00 # choose the new ratio for the comparison here

PN_area_ref = np.zeros((len(lat),len(lon)))
PN_bool_ref = np.where(bio_PN_tot[:,:] > ref_PN, 1, 0)
PN_A_ref = np.nansum(PN_bool_ref[:,:]*area,axis=(0,1))

PN_area_new = np.zeros_like(PN_area_ref)
PN_bool_new = np.where(bio_PN_tot[:,:] > new_PN, 1, 0)
PN_A_new = np.nansum(PN_bool_new[:,:]*area,axis=(0,1))

change_PN = (PN_A_ref-PN_A_new)/PN_A_ref
print(change_PN)

#%% Calculate differences in area for Fe:N of different values
new_FeN = 1.0

FeN_area_ref = np.zeros((len(lat),len(lon)))
FeN_bool_ref = np.where(bio_FeN_tot[:,:] > ref_FeN, 1, 0)
FeN_A_ref = np.nansum(FeN_bool_ref[:,:]*area,axis=(0,1))

FeN_area_new = np.zeros_like(FeN_area_ref)
FeN_bool_new = np.where(bio_FeN_tot[:,:] > new_FeN, 1, 0)
FeN_A_new = np.nansum(PN_bool_new[:,:]*area,axis=(0,1))

change_FeN = (FeN_A_ref-FeN_A_new)/FeN_A_ref
print(change_FeN)

#%% create masks for the altered P:N and Fe:N ratios to calculate the accuracy of correctly
# predicted points
mask_PN = np.where((bio_FeN_tot[:,:] > ref_FeN) & (bio_PN_tot[:,:] > new_PN), 1, 0) #keep Fe:N constant while varying P:N
mask_FeN = np.where((bio_FeN_tot[:,:] > new_FeN) & (bio_PN_tot[:,:] > ref_PN), 1, 0) #keep P:N constant while varying Fe:N

#%% calculate accuracy of prediction
list_idx = 1 #to chose which data from diaz_data_list to plot
# gives fraction of abundances that are within the predicted province
IN_PN = np.sum(mask_PN[lat_corr,lon_corr])/len(lat_corr)
IN_FeN = np.sum(mask_FeN[lat_corr,lon_corr])/len(lat_corr)
print(IN_PN)
print(IN_FeN)

acc_change_PN = (IN-IN_PN)/IN_PN
acc_change_FeN = (IN-IN_FeN)/IN_FeN
print(acc_change_PN)
print(acc_change_FeN)

#%% Create loop to vary P:N and Fe:N values and extract area and accuracy
# idea: plot P:N and Fe:N vs. area and accuracy

#P:N
new_PN = np.arange(0.5,1.5,0.01) # choose the new ratio for the comparison here

area_PN_var = np.zeros(len(new_PN))
change_PN_var = np.zeros_like(area_PN_var)
IN_PN_var = np.zeros_like(area_PN_var)
acc_PN_var = np.zeros_like(area_PN_var)

for i in range(len(new_PN)):
    PN_bool_new = np.where(bio_PN_tot[:,:] > new_PN[i], 1, 0)
    area_PN_var[i] = np.nansum(PN_bool_new[:,:]*area,axis=(0,1))
    change_PN_var[i] = (PN_A_ref-area_PN_var[i])/PN_A_ref
    #print(change_PN)
    mask_PN = np.where((bio_FeN_tot[:,:] > ref_FeN) & (bio_PN_tot[:,:] > new_PN[i]), 1, 0)
    IN_PN_var[i] = np.sum(mask_PN[lat_corr,lon_corr])/len(lat_corr)
    #print(IN_PN_var[i])
    acc_PN_var[i] = (IN-IN_PN_var[i])/IN_PN_var[i]
    #print(acc_PN_var[i])
    
#%% Create loop to vary P:N and Fe:N values and extract area and accuracy
# idea: plot P:N and Fe:N vs. area and accuracy

#Fe:N
new_FeN = np.arange(0.5,1.5,0.01) # choose the new ratio for the comparison here

area_FeN_var = np.zeros(len(new_FeN))
change_FeN_var = np.zeros_like(area_FeN_var)
IN_FeN_var = np.zeros_like(area_FeN_var)
acc_FeN_var = np.zeros_like(area_FeN_var)

for i in range(len(new_FeN)):
    FeN_bool_new = np.where(bio_FeN_tot[:,:] > new_FeN[i], 1, 0)
    area_FeN_var[i] = np.nansum(FeN_bool_new[:,:]*area,axis=(0,1))
    change_FeN_var[i] = (FeN_A_ref-area_FeN_var[i])/FeN_A_ref
    #print(change_FeN)
    mask_FeN = np.where((bio_FeN_tot[:,:] > ref_FeN) & (bio_FeN_tot[:,:] > new_FeN[i]), 1, 0)
    IN_FeN_var[i] = np.sum(mask_FeN[lat_corr,lon_corr])/len(lat_corr)
    #print(IN_FeN_var[i])
    acc_FeN_var[i] = (IN-IN_FeN_var[i])/IN_FeN_var[i]
    #print(acc_FeN_var[i])

#%% Vary both phi values at the same time
#Create loop to vary P:N and Fe:N values and extract area and accuracy
# idea: plot P:N and Fe:N vs. area and accuracy


#Fe:N
new_FeN = np.arange(0.5,1.5,0.01) # choose the new ratio for the comparison here
new_PN = np.arange(0.5,1.5,0.01)

area_var = np.zeros((len(new_FeN),len(new_PN)))
IN_var = np.zeros_like(area_var)

for i in range(len(new_FeN)):
    for j in range(len(new_PN)):
        #bool_new[i,j] = np.where(bio_FeN_tot[:,:] > new_FeN[i], 1, 0)
        mask_both = np.where((bio_FeN_tot[:,:] > new_FeN[i]) & (bio_PN_tot[:,:] > new_PN[j]), 1, 0)
        area_var[i,j] = np.nansum(mask_both[:,:]*area,axis=(0,1))
        IN_var[i,j] = np.sum(mask_both[lat_corr,lon_corr])/len(lat_corr)

#%%
fig,ax = plt.subplots(figsize=(12,6))
ax.imshow(new_FeN, new_PN, IN_var)

#%% Plot the differences in area and phi values
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

fig = plt.figure(figsize=(12,8))
ax = fig.gca(projection='3d')
# Plot the surface.
surf = ax.plot_surface(new_FeN, new_PN, IN_var[:,:])#, cmap=cm.coolwarm,linewidth=0, antialiased=False)

def randrange(n, vmin, vmax):
    '''
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    '''
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
    xs = randrange(n, 23, 32)
    ys = randrange(n, 0, 100)
    zs = randrange(n, zlow, zhigh)
    ax.scatter(xs, ys, zs, c=c, marker=m)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

#%% plot the results for this variation of supply ratios
fig,ax = plt.subplots(2,2,sharex=True,figsize=(9,6))
ax[0,0].plot(new_FeN[:],area_FeN_var[:],color='C0',linewidth=1.5)#,label='total area')
ax[1,0].plot(new_FeN[:],IN_FeN_var[:],color='C1',linewidth=1.5)#,label='22 mmol m$^{-3}$') 
ax[0,1].plot(new_PN[:],area_PN_var[:],color='C1',linewidth=1.5)#,label='22 mmol m$^{-3}$') 
ax[1,1].plot(new_PN[:],IN_PN_var[:],color='C2',linewidth=1.5)#,label='5 mmol m$^{-3}$') 
ax[1,0].set_xlabel('Fe:N',labelpad=10)
ax[1,1].set_xlabel('P:N',labelpad=10)
ax[0,0].set_ylabel('Area',labelpad=10)
ax[1,0].set_ylabel('Accuracy',labelpad=10)
#ax[0].set_yticks([0,4])
#ax[1].set_yticks([0,1.5])
#ax[2].set_yticks([0,0.2])
#ax[0].set_ylim(0,8)
#ax[1].set_ylim(0,3)
#ax[2].set_ylim(0,0.4)
#plt.suptitle('a)',x=0.08,y=0.95)
#plt.subplots_adjust(hspace=0)
plt.show()

#%% plot the results for this variation of supply ratios
fig,ax = plt.subplots(1,2,sharex=True,figsize=(9,3))
p1 = ax[0].plot(new_FeN[:],area_FeN_var[:],color='C0',linewidth=1.5,label='varying Fe:N')
p2 = ax[1].plot(new_FeN[:],IN_FeN_var[:],color='C0',linewidth=1.5,label='varying Fe:N') 
ax[0].plot(new_PN[:],area_PN_var[:],color='C1',linewidth=1.5,label='varying P:N')
ax[1].plot(new_PN[:],IN_PN_var[:],color='C1',linewidth=1.5,label='varying P:N')
ax[0].legend(loc='lower left')
ax[1].legend(loc='lower left')
ax[0].set_xlabel('Fe:N, P:N')
ax[0].set_ylabel('Area')
ax[1].set_ylabel('Accuracy')

plt.show()

# thoughts on this:
# why don't we get an accuracy of 100% with very small phi? --> because the total area
# is calculated for Fe:N & P:N. Thus, at some point the area does not get bigger and the
# accuracy not higher because the other nutrient is limiting (the one that's kept constant)
# --> vary both nutrients at the same time

#%% Calculate absences --> meaning obs - nz_obs
#absences = np.where(find_obs[0]==1 & find_nz_obs[0]==0)
#absence = np.where(find_obs[0][:] != find_nz_obs[0][:])

#maybe write a loop?
#absence = np.zeros_like(obs_tot)
#for i in range(len(find_obs[0])):
#    if find_obs[0][i] == find_nz_obs[0][i]:
#        absence[0,i] = 0
#    else:
#        absence[0,i] = 1
    
#%% just a plot to quickly display variables

col = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(12,4))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,mask_PN,levels=np.linspace(-1,2,20),cmap=col)#,extend='both')
#plt.plot(listi[1],listi[0],'.',color='m')
#plt.plot(pz_d[0],pz_d[1],'.',color='g')
#plt.plot(lon_d[find_obs[1]],lat_d[find_obs[0]],'.',color='b')
plt.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='g')
#plt.plot(lon_d[find_bm[1]],lat_d[find_bm[0]],'.',color='r')
#plt.plot(lon_d[find_nifH[1]],lat_d[find_nifH[0]],'.',color='c')
#plt.plot(lon_d[find_nz_obs[1]],lat_d[find_nz_obs[0]],'.',color='m')
#plt.plot(lon_d[find_nz_abund[1]],lat_d[find_nz_abund[0]],'.',color='orange')
#plt.plot(lon_d[find_nz_bm[1]],lat_d[find_nz_bm[0]],'.',color='k')
#plt.plot(lon_d[find_nz_nifH[1]],lat_d[find_nz_nifH[0]],'.',color='w')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
#cbar = plt.colorbar(c,ax=ax)
#cbar.set_label('transport '+str(label_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()

#%% Plot 1 nutrient at 1 depth

nu = 1   #chose nutrient here: 0=Fe, 1=N, 2=P
nutr = [bio_PN_tot,bio_FeN_tot,bio_PN_trans,bio_FeN_trans,bio_PN_remin,bio_FeN_remin]
label_nut = ['N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)','Fe (mol/m$^{2}$/y)']
name_nut = ['P:N','Fe:N','transport P:N','transport Fe:N','remin P:N','remin Fe:N']

levs_PN = np.linspace(0,1000,11)
levs_FeN = np.linspace(-2,2,21)
levs = [levs_PN,levs_FeN,levs_PN,levs_FeN,levs_PN,levs_FeN]

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(12,4))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,np.log(nutr[nu]),levels=levs[nu],cmap=colmap,extend='both')
con1 = ax.contour(lon,lat,nutr[nu],levels=[1.2],colors='k',linewidths=1,linstyle='solid')
con2 = ax.contour(lon,lat,nutr[nu],levels=[2.5],colors='r',linewidths=1,linstyle='solid')

#plt.plot(lon_d[find_obs[1]],lat_d[find_obs[0]],'.',color='b')
#plt.plot(lon_d[find_abund[1]],lat_d[find_abund[0]],'.',color='g')
#plt.plot(lon_d[find_bm[1]],lat_d[find_bm[0]],'.',color='r')
#plt.plot(lon_d[find_nifH[1]],lat_d[find_nifH[0]],'.',color='c')
#plt.plot(lon_d[find_nz_obs[1]],lat_d[find_nz_obs[0]],'.',color='m')
#plt.plot(lon_d[find_nz_abund[1]],lat_d[find_nz_abund[0]],'.',color='orange')
#plt.plot(lon_d[find_nz_bm[1]],lat_d[find_nz_bm[0]],'.',color='k')
#plt.plot(lon_d[find_nz_nifH[1]],lat_d[find_nz_nifH[0]],'.',color='w')

lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label(''+str(name_nut[nu])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/darfts_trans_NO3_fixed.png', bbox_inches='tight', dpi=300)


#%% Plot 1 nutrient at 1 depth

nut = 0    #chose nutrient here: 0=Fe, 1=N, 2=P
nutrient = [bio_PN_tot,bio_PN_trans,bio_PN_remin]
label_nut = ['N (mol/m$^{2}$/y)','P (mol/m$^{2}$/y)','Fe (mol/m$^{2}$/y)']
name_nut = ['P:N','transport P:N','remin P:N']

levs_PN = np.linspace(0.5,1.5,21)
levs_FeN = np.linspace(-2,2,21)
levs = [levs_PN,levs_FeN]

colmap = plt.get_cmap('RdBu_r')

fig,ax = plt.subplots(subplot_kw={'projection':ccrs.PlateCarree(central_longitude=0)},figsize=(12,4))
ax.coastlines(color='#888888',linewidth=1.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m', edgecolor='none', facecolor=cfeature.COLORS['land']))
c = ax.contourf(lon,lat,nutrient[nut],levels=levs[nut],cmap=colmap,extend='both')
#con1 = ax.contour(lon,lat,nutrient[nut],levels=[0.99],colors='purple',linewidths=1,linstyle='solid')
con2 = ax.contour(lon,lat,nutrient[nut],levels=[1.04],colors='r',linewidths=1,linstyle='solid')
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks([0,60,120,180,240,300,360], crs=ccrs.PlateCarree())
ax.set_yticks([-80, -60, -40, -20, 0, 20, 40, 60, 80], crs=ccrs.PlateCarree())
cbar = plt.colorbar(c,ax=ax)
cbar.set_label(''+str(name_nut[nut])+'',rotation=90, position=(0.5,0.5))
plt.show()
#fig.savefig('/Users/meilers/MITinternship/Plots/overview_nutr_bioav_'+str(name_nut[nut])+'_104.png', bbox_inches='tight', dpi=300)

#%% Calculate area elementwise --> yields same result as from section above

#PN_A_r = np.zeros((160,360))
#for i in range(len(lat)):
#    for j in range(len(lon)):
#        PN_A_r[i,j] = PN_bool_ref[i,j]*area[i,j]
#PN_tot_r = np.sum(PN_A_r)    
#
#PN_A_n = np.zeros((160,360))
#for i in range(len(lat)):
#    for j in range(len(lon)):
#        PN_A_n[i,j] = PN_bool_new[i,j]*area[i,j]
#PN_tot_n = np.sum(PN_A_n)
#
#change2 = (PN_tot_r-PN_tot_n)/PN_tot_r
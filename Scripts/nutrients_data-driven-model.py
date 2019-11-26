# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cmocean as cm

#%% Load data

# iron dust flux
fein = np.fromfile('/Users/meilers/MITinternship/Data/mahowald2009_solubile_current_smooth_oce_mth-2d.bin', '>f4').reshape(12,160,360)

grid = xr.open_dataset('/Users/meilers/MITinternship/Data/supply50m.nc')
area_info = xr.open_dataset('/Users/meilers/MITinternship/Data/grid.nc')

lon = grid.lon   #needed for plotting
lat = grid.lat    #needed for plotting

area = area_info.rA

#%% Read in observations of Tang and Cassar, 2019
diazotroph_observations = pd.read_csv(r'/Users/meilers/MITinternship/Data/Tang_and_Cassar-2019/nifH_Gene_Integral.csv')
#print(diazotroph_observations)

longitude = pd.DataFrame(diazotroph_observations, columns = ['LONGITUDE'])
latitude = pd.DataFrame(diazotroph_observations, columns = ['LATITUDE'])

nifH_Tri = pd.DataFrame(diazotroph_observations, columns = ['Trichodesmium nifH Gene (x106 copies m-2)'])
nifH_Tri_bm = pd.DataFrame(diazotroph_observations, columns = ['Trichodesmium Biomass (mg C m-2)'])
nifH_UCYN_A = pd.DataFrame(diazotroph_observations, columns = ['UCYN-A nifH Gene (x106 copies m-2)'])
nifH_UCYN_A_bm = pd.DataFrame(diazotroph_observations, columns = ['UCYN-A Biomass Conversion factor (mg C/106 nifH copies)'])

#nifH_Tri = diazotroph_observations['Trichodesmium nifH Gene (x106 copies m-2)'].values.tolist()

#%%
plt.scatter(longitude, latitude, c=nifH_Tri[0:10],cmap='viridis')
plt.show()

#%%
ax2 = diazotroph_observations.plot.scatter(x='LONGITUDE',
                                           y='LATITUDE',
                                           c=nifH_Tri,
                                           colormap='viridis')

#%% 
mask = np.zeros_like(DIN)
mask[nifH_Tri>0] = 1

#%% Read in monthly nutrient data

months_vec = range(0,12)
#months_list = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mon_list = ['26160','26400','26640','26880','27120','27360','27600','27840','28080','28320','28560','28800']

#timevector = np.zeros((1,1))
# create empty arrays to next fill with the nutrient data from the model
NH4 = np.zeros((len(months_vec),6,160,360))
NO2 = np.zeros((len(months_vec),6,160,360))
NO3 = np.zeros((len(months_vec),6,160,360))
DIP = np.zeros((len(months_vec),6,160,360))

# Open dataset and load values for each nutrient (top 100m --> 0:6 of depth dimension)
i = 0
for month in months_vec:
    data = xr.open_dataset('/Users/meilers/MITinternship/Data/run19_33_MONTH/3d.00000'+str(mon_list[month])+'.nc')
    NH4[month,:,:,:] = data.TRAC02.values[:,0:6,:,:]
    NO2[month,:,:,:] = data.TRAC03.values[:,0:6,:,:]
    NO3[month,:,:,:] = data.TRAC04.values[:,0:6,:,:]
    DIP[month,:,:,:] = data.TRAC05.values[:,0:6,:,:]
    data.close()
    print('Month: '+str(month)) #as control if loop works
    i += 1

DIN = NH4 + NO2 + NO3


#%% If only the surface nutrient concentration is of interest
months_vec = range(0,12)
#months_list = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
mon_list = ['26160','26400','26640','26880','27120','27360','27600','27840','28080','28320','28560','28800']

#timevector = np.zeros((1,1))
# create empty arrays to next fill with the nutrient data from the model
NH4 = np.zeros((len(months_vec),160,360))
NO2 = np.zeros((len(months_vec),160,360))
NO3 = np.zeros((len(months_vec),160,360))
DIP = np.zeros((len(months_vec),160,360))

# Open dataset and load values for each nutrient (top 100m --> 0:6 of depth dimension)
i = 0
for month in months_vec:
    data = xr.open_dataset('/Users/meilers/MITinternship/Data/run19_33_MONTH/3d.00000'+str(mon_list[month])+'.nc')
    NH4[month,:,:] = data.TRAC02.values[:,0,:,:]
    NO2[month,:,:] = data.TRAC03.values[:,0,:,:]
    NO3[month,:,:] = data.TRAC04.values[:,0,:,:]
    DIP[month,:,:] = data.TRAC05.values[:,0,:,:]
    data.close()
    print('Month: '+str(month)) #as control if loop works
    i += 1

DIN = NH4 + NO2 + NO3

DIN_int = np.mean(DIN,axis=0) 
DIP_int = np.mean(DIP,axis=0)

#%% Read in diazotroph data from model - same setup as for nutrients
    # 5 diazotroph species. according to Steph TRAC30 to TRAC34

diaz1 = np.zeros((len(months_vec),6,160,360))
diaz2 = np.zeros((len(months_vec),6,160,360))
diaz3 = np.zeros((len(months_vec),6,160,360))
diaz4 = np.zeros((len(months_vec),6,160,360))
diaz5 = np.zeros((len(months_vec),6,160,360))

# Open dataset
i = 0
for month in months_vec:
    diaz = xr.open_dataset('/Users/meilers/MITinternship/Data/run19_33_MONTH/3d.00000'+str(mon_list[month])+'.nc')
    diaz1[month,:,:,:] = diaz.TRAC30.values[:,0:6,:,:]
    diaz2[month,:,:,:] = diaz.TRAC31.values[:,0:6,:,:]
    diaz3[month,:,:,:] = diaz.TRAC32.values[:,0:6,:,:]
    diaz4[month,:,:,:] = diaz.TRAC33.values[:,0:6,:,:]
    diaz5[month,:,:,:] = diaz.TRAC34.values[:,0:6,:,:]
    diaz.close()
    print('Month: '+str(month))
    i += 1
# Sum up diazotroph data into one array
diaz = diaz1 + diaz2 + diaz3 + diaz4 + diaz5

#%% define sec per year and dz
sec = 1  #seconds per year (365.25*86400) - not needed; I used it for the annual analysis. 
dz = [10,10,15,20,20,25] #,35,50,75,100] #(...) dz between two depth layers
#dz = [10,0,0,0,0,0]
#%% sum nutrients up over depth and multiply with corresponding dz and sec

DIN_int = np.zeros((12,6,160,360))              #create empty vector
for i in range(len(dz)):                        #loop over depths                 
    DIN_int[:,i,:,:] = DIN[:,i,:,:]*dz[i]*sec   #multiply fluxes at each depth with corresponding dz
    print(np.max(DIN_int[:,i,:,:]))
    #print(i)
    #print(dz[i])
DIN_int = np.mean(DIN_int,axis=(0,1))                #sum up over depth --> new shape of array (12,160,360)

DIP_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    DIP_int[:,i,:,:] = DIP[:,i,:,:]*dz[i]*sec    
DIP_int = np.mean(DIP_int,axis=(0,1))

diaz_int = np.zeros((12,6,160,360))
for i in range(len(dz)):
    diaz_int[:,i,:,:] = diaz[:,i,:,:]*dz[i]*sec    
diaz_int = np.sum(diaz_int,axis=1)
diaz_int = np.mean(diaz_int,axis=0)*10 #multiply with 10 to convert from 100m-2 to l-1; no conversion of biomass to nifH genes yet.
  
#%% Calculate N:P

N_P = np.log(DIN_int)/np.log(DIP_int)

#%% Plot relations between environmental properties (DIN, DIP, N:P) and volumetric diazotroph abundances

fig,ax = plt.subplots(1,1,figsize=(4,4))
plt.plot(np.log(DIN_int),np.log(diaz_int),'.')
plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 17:11:35 2020

@author: meilers
"""
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from skimage import measure

#%% Create new mask of regions from Darwin model output

# Read in monthly nutrient data

months_vec = range(0,12)
mon_list = ['26160','26400','26640','26880','27120','27360','27600','27840','28080','28320','28560','28800']
sec = 31557600

# create empty arrays to next fill with the nutrient data from the model
PO4 = np.zeros((len(months_vec),160,360))
S_PO4 = np.zeros((len(months_vec),160,360))

# Open dataset and load values for each nutrient (top 100m --> 0:6 of depth dimension)
i = 0
for month in months_vec:
    data = xr.open_dataset('/Users/meilers/MITinternship/Data/run19_33_MONTH/nutr_tend.00000'+str(mon_list[month])+'.nc')
    PO4[month,:,:] = data.gTr05.values[:,0,:,:]
    S_PO4[month,:,:] = data.S_PO4.values[:,0,:,:]
    data.close()
    print('Month: '+str(month)) #as control if loop works
    i += 1

PO4_t = np.add(PO4,S_PO4)   # sum the two PO4 sources up 
PO4_tot = np.mean(PO4_t,axis=0)*sec # create annual mean of total PO4

thresh = 0.3 #mmol/m3

#PO4_tot_above = np.where((PO4_tot > thresh), 1, 0)
#PO4_tot_below = np.where((PO4_tot < thresh), 1, 0)

#%% create the regions with the measure.label function from skimage

### HOW CAN I SET THE CONTINENTS TO BACKGROUND?
### HOW CAN I INCORPORATE THE THRESHOLD?

mask_darwin = measure.label(PO4_tot, neighbors=None, background=None, return_num=False, connectivity=None)

#%% Just for some quick plots
fig,ax = plt.subplots(figsize=(6,4))
c = ax.contourf(mask_darwin)
cbar = plt.colorbar(c,ax=ax)
plt.show()
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 10:25:50 2024

@author: pb130
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt



dt = 7*24*3600 
# nc = Dataset('present_day_without_ozone_hole.nc','r')
# nc = Dataset('Holocene_simulation/Holocene_base_case_1980_oz.nc','r')
# nc = Dataset('Glacial_simulation/Glacial_base_case_1980_oz.nc','r')
nc = Dataset('Glacial_simulation/Glacial_base_case_1980_oz_with_phiT.nc','r')

dt = 7*24*3600 

fout = nc['parameters']['fout'][:]
weekly_profile = nc['outputs2']
# weekly variables
FD_mass = weekly_profile['FD_mass'][:]*dt
FP_mass = weekly_profile['FP_mass'][:]*dt
depth_profile = nc['outputs1']
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_data = depth_profile['SN_data'][:]
SN_thickness = depth_profile['SN_thickness'][:]
# nc.close()

# calculate the remaining fraction from FA and FD
nb_tt = FD_mass.shape[0]
FA_mass = np.zeros((nb_tt))
FA_d15N = np.zeros((nb_tt))
FA_D17O = np.zeros((nb_tt))

SN_d = 350 # snow density
Fpri = 2.2e-6

for i in range(nb_tt):
    idx = np.where(SN_data==i) # find the index
    FA = SN_thickness*SN_d*SN_conc*1e-9*14/62/dt # data should be date here. Ooop.
    FA_mass[i] = FA[idx].sum()
    FA_d15N[i] = (FA[idx]*SN_d15N[idx]).sum()/FA[idx].sum()
# using the annual value in 5th year
# idx = 260: 312
FA_annual = FA_mass[260: 312].mean()*52*dt
FA_d15N_annual = (FA_mass[260: 312]*FA_d15N[260: 312]).sum()/FA_mass[260: 312].sum()

print('%.2f  %.2f'%(FA_d15N_annual,(1-FA_annual/Fpri)*100))
    





























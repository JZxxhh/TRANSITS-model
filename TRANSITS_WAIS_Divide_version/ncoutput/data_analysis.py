# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:40:23 2020

@author: admin
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


SN_ob = np.genfromtxt('WD_snowpack_nitrate.txt',skip_header=1)



dt = 7*24*3600 
# nc = Dataset('test_SE_Dome.nc','r')
nc = Dataset('present_day_with_ozone_hole.nc','r')

fout = nc['parameters']['fout'][:]
weekly_profile = nc['outputs2']

# weekly variables
FP_mass1 = weekly_profile['FP_mass'][:]*dt
FP_d15N = weekly_profile['FP_d15N'][:]
FP_D17O = weekly_profile['FP_D17O'][:]
FD_mass = weekly_profile['FD_mass'][:]*dt
FD_d15N = weekly_profile['FD_d15N'][:]
FD_D17O = weekly_profile['FD_D17O'][:]
eps_surface = weekly_profile['eps_surface'][:]
eps_surface = weekly_profile['eps_surface'][:]
frac_rem1 = weekly_profile['frac_rem'][:][:-1]
surface_conc = weekly_profile['surface_conc'][:-25]
surface_d15N = weekly_profile['surface_d15N'][:-25]
surface_D17O = weekly_profile['surface_D17O'][:-25]

depth_profile = nc['outputs1']
SN_conc1 = depth_profile['SN_conc'][:]
SN_d15N1 = depth_profile['SN_d15N'][:]
SN_data = depth_profile['SN_data'][:]
SN_D17O = depth_profile['SN_D17O'][:]
SN_thickness = depth_profile['SN_thickness'][:]
nc.close()

depth = SN_thickness.cumsum()


nc = Dataset('present_day_without_ozone_hole.nc','r')

fout = nc['parameters']['fout'][:]
weekly_profile = nc['outputs2']
FP_mass2 = weekly_profile['FP_mass'][:]*dt
frac_rem2 = weekly_profile['frac_rem'][:][:-1]

depth_profile = nc['outputs1']
SN_conc2 = depth_profile['SN_conc'][:]
SN_d15N2 = depth_profile['SN_d15N'][:]
nc.close()



# plt.plot(SN_ob[:,0]/100,SN_ob[:,1])


plt.plot(depth, SN_d15N1)
plt.plot(depth, SN_d15N2)
# plt.plot(depth,frac_rem2)

plt.xlim([0,2])

# calculate 2-m avergae d15N in comparison with observation
# idx = np.searchsorted(depth, 2)
# d15N_2m =  (SN_conc1[:idx]*SN_d15N1[:idx]).sum()/SN_conc1[:idx].sum()







# # calculate the remaining fraction from FA and FD
# nb_tt = FD_mass.shape[0]
# FA_mass = np.zeros((nb_tt))
# FA_d15N = np.zeros((nb_tt))
# FA_D17O = np.zeros((nb_tt))

# SN_d = 350 # snow density
# Fpri = 2.7e-6

# for i in range(nb_tt):
#     idx = np.where(SN_data==i) # find the index
#     FA = SN_thickness*SN_d*SN_conc*1e-9*14/62/dt # data should be date here. Ooop.
#     FA_mass[i] = FA[idx].sum()
#     FA_d15N[i] = (FA[idx]*SN_d15N[idx]).sum()/FA[idx].sum()
#     FA_D17O[i] = (FA[idx]*SN_D17O[idx]).sum()/FA[idx].sum() # note D17O is inaccurate here becasue I didn't modify the oxidants in model.

# # using the annual value in 5th year

# # idx = 260: 312
# FA_annual = FA_mass[260: 312].mean()*52*dt
# FA_d15N_annual = (FA_mass[260: 312]*FA_d15N[260: 312]).sum()/FA_mass[260: 312].sum()

# print('%.2f  %.2f'%(FA_d15N_annual,(1-FA_annual/Fpri)*100))














































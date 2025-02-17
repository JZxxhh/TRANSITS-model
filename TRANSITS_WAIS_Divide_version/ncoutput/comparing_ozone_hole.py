# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 20:23:10 2024

@author: pb130
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


dt = 7*24*3600 
nc = Dataset('present_day_without_ozone_hole.nc','r')

depth_profile = nc['outputs1']
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_thickness = depth_profile['SN_thickness'][:]
SN_data = depth_profile['SN_data'][:]
weekly_profile = nc['outputs2']
FD_mass = weekly_profile['FD_mass'][:]
nc.close()

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

SN_d15N_without_oz = SN_d15N
SN_d15N_without_oz_ave = (FA_mass[260: 312]*FA_d15N[260: 312]).sum()/FA_mass[260: 312].sum() # archived average

nc = Dataset('present_day_with_ozone_hole.nc','r')

depth_profile = nc['outputs1']
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_thickness = depth_profile['SN_thickness'][:]
SN_data = depth_profile['SN_data'][:]
weekly_profile = nc['outputs2']
FD_mass = weekly_profile['FD_mass'][:]
nc.close()

nb_tt = FD_mass.shape[0]
FA_mass = np.zeros((nb_tt))
FA_d15N = np.zeros((nb_tt))
FA_D17O = np.zeros((nb_tt))

for i in range(nb_tt):
    idx = np.where(SN_data==i) # find the index
    FA = SN_thickness*SN_d*SN_conc*1e-9*14/62/dt # data should be date here. Ooop.
    FA_mass[i] = FA[idx].sum()
    FA_d15N[i] = (FA[idx]*SN_d15N[idx]).sum()/FA[idx].sum()

SN_d15N_with_oz = SN_d15N
SN_d15N_with_oz_ave = (FA_mass[260: 312]*FA_d15N[260: 312]).sum()/FA_mass[260: 312].sum() # archived average
depth = SN_thickness.cumsum()

gs_kw = dict(height_ratios=[2, 1])
fig, axd = plt.subplot_mosaic([['upper'],
                               ['lower']],sharex=True,
                              gridspec_kw=gs_kw, figsize=(6,4),
                              layout="constrained")

plt.subplots_adjust(hspace=0)

axd['upper'].plot(depth,SN_d15N_with_oz,'y-',
                  label=r'With ozone hole')
axd['upper'].plot(depth,SN_d15N_without_oz,'b-',
                  label=r'Without ozone hole')
axd['upper'].set_ylabel(r'$\delta^{15}$N(NO$_{3}$$^{-}$)',fontsize=14)
axd['upper'].set_xlim([0,2])

axd['upper'].legend(loc='lower center',fontsize=10)

difference = SN_d15N_with_oz-SN_d15N_without_oz
axd['lower'].plot(depth,difference,'k-',
                  label=r'$\delta^{15}$N(NO$_{3}$$^{-}$) difference')

axd['lower'].set_xlabel('Depth (m)',fontsize=14)

axd['lower'].set_ylabel(r'Difference',fontsize=14)

axd['upper'].text(0.02,0.9,'(a)',fontsize=14,transform=axd['upper'].transAxes)
axd['lower'].text(0.02,0.8,'(b)',fontsize=14,transform=axd['lower'].transAxes)


plt.savefig('Comparing_d15N_oz_hole.png',dpi=300,bbox_inches='tight')







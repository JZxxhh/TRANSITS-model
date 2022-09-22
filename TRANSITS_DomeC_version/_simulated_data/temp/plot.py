# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 18:58:33 2022

@author: lenovo
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


file_path = 'TRANSITS_simul___.nc'
nc = Dataset(file_path,'r')
gp = nc.groups['outputs']

# snowpack depth profile
depth = gp['SN_depth_mid'][:]
conc = gp['SN_conc'][:]
d15N = gp['SN_d15N'][:]
D17O = gp['SN_D17O'][:]

nc.close()

# This file is used to illustrate the numerical diffusion casued by the
# RESAMPLE function in TRANSITS.
# There is clear diffusion of snow nitrate concentration and isotops casued 
# even without turning on the diffusion process.
# This effect can be very significant if the accumulated snow thickness in each
# timestep satisfies some specified conditions.

###### PLOT ######
fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(9, 5))
plt.subplots_adjust(hspace=0.1)
plt.subplots_adjust(wspace=0.1)

ax1.set_xticks([])
ax11 = ax1.twiny()
ax11.plot(conc[:,26], -depth[:,26], color='k',linestyle='-',linewidth=1,
         markersize=4,label='model simulation')

# ax11.set_ylim([-0.5,0])

ax2.set_xticks([])
ax21 = ax2.twiny()
ax21.plot(d15N[:,26], -depth[:,26], color='k',linestyle='-',linewidth=1,
         markersize=4,label='model simulation')
# ax21.set_ylim([-0.5,0])

ax3.set_xticks([])
ax31 = ax3.twiny()
ax31.plot(D17O[:,26], -depth[:,26], color='k',linestyle='-',linewidth=1,
         markersize=4,label='model simulation')
# ax31.set_ylim([-0.5,0])

ax11.set_xlabel(r'[NO$_{3}$$^{-}$], ng.g$^{-1}$', fontsize=14, color='k', ha="center", va="center", labelpad=14)
ax21.set_xlabel(r'10$^{3} \times$ $\delta^{15}$N', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel
ax31.set_xlabel(r'10$^{3} \times$ $\Delta^{17}$O', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel
ax1.set_ylabel('depth (m)', fontsize=14, color='k', ha="center", va="center", labelpad=14)

plt.savefig('DomeC_result.png',dpi=300,bbox_inches='tight')















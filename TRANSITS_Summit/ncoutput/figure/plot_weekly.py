# -*- coding: utf-8 -*-
"""
Created on Wed Dec 04 05:26:08 2019

@author: admin
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

dt = 7*24*3600 
nc = Dataset('Summit_PD_fexp_0.35.nc','r')
depth_profile = nc['outputs1']
weekly_profile = nc['outputs2']

# weekly variables
FP_mass = weekly_profile['FP_mass'][:156]*dt
FP_d15N = weekly_profile['FP_d15N'][:156]
FP_D17O = weekly_profile['FP_D17O'][:156]
FD_mass = weekly_profile['FD_mass'][:156]*dt
FD_d15N = weekly_profile['FD_d15N'][:156]
FD_D17O = weekly_profile['FD_D17O'][:156]
eps_surface = weekly_profile['eps_surface'][:156]
eps_surface = weekly_profile['eps_surface'][:156]
frac_rem = weekly_profile['frac_rem'][:156]
nc.close()


fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True,figsize=(6,9))
plt.subplots_adjust(wspace =0, hspace =0)

ax1.plot(np.arange(156),FP_mass,'r-', lw=2, label='FP')
ax1.plot(np.arange(156),FD_mass,'k-', lw=2, label='FD')
ax1.text(10,1.92e-7,'(a)',fontsize=14)
ax1.text(-10,2.15e-7,'10$^{-7}$',fontsize=10)
ax1.set_yticks(np.arange(0,2.2,0.3)*1e-7)
#ax1.set_ylim([0,2.2e-7])
ax1.set_ylabel('kgN.m$^{-2}$.week$^{-1}$',fontsize=14,color='k')

axy = ax1.twiny()
axy.set_xticks([0,1,2,3])
axy.set_xticklabels('')
axy.set_xticks([0.5,1.5,2.5], minor=True)
axy.set_xticklabels([r'year 1$^{st}$',r'year 2$^{nd}$',r'year 3$^{rd}$'], minor=True,fontsize=14)
ax1.set_yticklabels([0,0.3,0.6,0.9,1.2,1.5,1.8,2.1])

ax11 = ax1.twinx()
ax11.plot(np.arange(156),FP_mass/FD_mass,linestyle='--',color='b')
ax1.legend(loc=(0.27,0.75))
ax11.set_ylabel('FP / FD',fontsize=14,color='b')


FP_d15N[np.where(FP_d15N==0)] = np.nan
ax2.plot(np.arange(156),FP_d15N,'r-', lw=2, label='FP')
ax2.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N',fontsize=14,color='r')
ax2.text(10,-66,'(b)',fontsize=14)
ax22 = ax2.twinx()
ax22.plot(np.arange(156),FD_d15N,'k-', lw=2, label='FD')
ax22.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N',fontsize=14,color='k')

FP_D17O[np.where(FP_D17O==0)] = np.nan
ax3.plot(np.arange(156),FP_D17O,'r-', lw=2, label='FP')
ax3.set_ylabel(r'10$^{3} \times$ $\Delta^{17}$O',fontsize=14,color='r')
ax33 = ax3.twinx()
ax33.plot(np.arange(156),FD_D17O,'k-', lw=2, label='FD')
ax33.set_ylabel(r'10$^{3} \times$ $\Delta^{17}$O',fontsize=14,color='k')
ax3.set_xlabel('week of year',fontsize=14)
#ax3.set_xlim([0,156])
ax33.xaxis.set_major_locator(ticker.FixedLocator((np.arange(7)*26)))
ax33.xaxis.set_major_formatter(ticker.FixedFormatter(([0,26,0,26,0,26,0])))

ax3.text(10,33,'(c)',fontsize=14)

plt.savefig('weekly_flux.png', bbox_inches='tight', dpi=300)
# plt.savefig("Figure 2.pdf", bbox_inches='tight')













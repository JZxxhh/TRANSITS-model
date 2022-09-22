# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 21:25:48 2020

@author: admin
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def cal_mean(interval,depth,conc,d15N,D17O,end=2.1):
    sp = int(end/interval)
    dp = np.linspace(0,end,sp+1)
    idx = np.searchsorted(depth,dp)
    res = np.zeros((len(dp)-1,4))
    for i in range(len(dp)-1):
        res[i,0] = dp[i:i+2].mean()
        res[i,1] = conc[idx[i]:idx[i+1]].mean()
        res[i,2] = (conc[idx[i]:idx[i+1]]*d15N[idx[i]:idx[i+1]]).sum()/conc[idx[i]:idx[i+1]].sum()
        res[i,3] = (conc[idx[i]:idx[i+1]]*D17O[idx[i]:idx[i+1]]).sum()/conc[idx[i]:idx[i+1]].sum()
    return res

# nc = Dataset('Summit_basecase_Javis.nc','r')
nc = Dataset('Summit_PD_fexp_0.35.nc','r')
# nc = Dataset('Summit_temperature_dependent_phi.nc','r')
depth_profile = nc['outputs1']
weekly_profile = nc['outputs2']

depth = depth_profile['SN_depth'][:]
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_D17O = depth_profile['SN_D17O'][:]
SN_data = depth_profile['SN_data'][:]
SN_thickness = depth_profile['SN_thickness'][:]

# weekly variables
FP_mass = weekly_profile['FP_mass'][:]
FP_d15N = weekly_profile['FP_d15N'][:]
FP_D17O = weekly_profile['FP_D17O'][:]
FD_mass = weekly_profile['FD_mass'][:]
FD_d15N = weekly_profile['FD_d15N'][:]
FD_D17O = weekly_profile['FD_D17O'][:]
eps_surface = weekly_profile['eps_surface'][:]
frac_rem = weekly_profile['frac_rem'][:]
nc.close()

#calculate the d15N of preserved snow nitrate
SN_d15N_week = np.zeros((len(FD_d15N)))
for i in range(len(FD_d15N)):
    idx1 = 0
    idx2 = 0
    for j in range(len(depth)):
        if SN_data[j] == i and SN_data[j-1] != i:
            idx1 = j
        if SN_data[j] != i and SN_data[j-1] == i:
            idx2 = j
            break
    SN_d15N_week[i] = (SN_conc[idx1:idx2]*SN_d15N[idx1:idx2]).sum()/SN_conc[idx1:idx2].sum()

PIF = SN_d15N_week-FD_d15N # PHOTOLYSIS INDUCED FRACTIONATION

fig,(ax,ax0) = plt.subplots(nrows=2,sharex=True,figsize=(8,8))
plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)

ax.set_ylim([-5,20])
ax.set_xlim([0,190])
# season range: spring: 9-21; autumn: 22-34; autumn: 35-47 winter: 48-8
# season index: 8, 20, 33, 47
# spring
for i in range(4):
    ax.fill([8+52*i,21+52*i,21+52*i,8+52*i],[-5,-5,20,20],'b',alpha=0.2)
# summer
for i in range(4):
    ax.fill([21+52*i,33+52*i,33+52*i,21+52*i],[-5,-5,20,20],'r',alpha=0.2)
# autumn
for i in range(3):
    ax.fill([33+52*i,47+52*i,47+52*i,33+52*i],[-5,-5,20,20],'g',alpha=0.2)
# plot average PIF
# index of season bin
idx = [0,9,22,35,48,61,74,87,100,113,126,139,152,165,178,187]
PIF_ave = np.zeros((len(idx)-1))
PIF_std = np.zeros((len(idx)-1))
for i in range(len(PIF_ave)):
    PIF_ave[i] = PIF[idx[i]:idx[i+1]].mean()
    PIF_std [i]= PIF[idx[i]:idx[i+1]].std()
x = (np.array(idx[1:])+(np.array(idx[:-1])))/2
# ax.errorbar(x,PIF_ave,yerr=PIF_std,fmt='o',ecolor='k',elinewidth=1,ms=5,mfc='wheat',mec='salmon',capsize=3,label='PIE_model') 
# data from Javis 2008
xob = [78,  91.5, 104.5, 117.5, 130.5, 143.5, 156.5, 169.5]
yob = [6.7, -3.13, 6.1, 10.22, 6.88, 1.46, 5.49, 18.37]
ax.scatter(xob,yob,s=26,c='k',marker='*',label='PIE_ob')

ax.xaxis.set_major_locator(ticker.FixedLocator((np.arange(8)*26)))
ax.xaxis.set_major_formatter(ticker.FixedFormatter(([0,26,0,26,0,26,0,26])))
ax.set_xlabel('week of year',fontsize=14)
ax.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N',fontsize=14)
ax.plot(np.arange(len(FD_d15N)),PIF,color='k',label='PIE')

ax1 = ax.twinx()
ax1.plot(np.arange(len(FD_d15N)),(1-frac_rem)*100,'r--',label=r'f$_{loss}$')
ax1.set_ylabel('loss fraction (%)',fontsize=14,color='r')

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax1.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2,loc='upper center',fontsize=10,ncol=4)
ax1.tick_params(axis='y', colors='r')
ax2 = ax.twiny()
ax2.set_xlim([0,190])
ax2.xaxis.set_major_locator(ticker.FixedLocator((np.arange(8)*26)))
ax2.xaxis.set_major_formatter(ticker.FixedFormatter((['','2004','','2005','','2006','','2007'])))
xticks = ax2.xaxis.get_major_ticks()
#for i in [1,3,5,7]:
#    xticks[i].label1.set_visible(False)
for i in range(4):
    ax0.fill([8+52*i,21+52*i,21+52*i,8+52*i],[-5,-5,20,20],'b',alpha=0.2)
# summer
for i in range(4):
    ax0.fill([21+52*i,33+52*i,33+52*i,21+52*i],[-5,-5,20,20],'r',alpha=0.2)
# autumn
for i in range(3):
    ax0.fill([33+52*i,47+52*i,47+52*i,33+52*i],[-5,-5,20,20],'g',alpha=0.2)
ax0.set_xlim([0,190])
ax0.plot(np.arange(len(FD_d15N)),PIF+FD_d15N,color='k',label='PIE')


# plt.savefig('Figure2_1.png',dpi=400)




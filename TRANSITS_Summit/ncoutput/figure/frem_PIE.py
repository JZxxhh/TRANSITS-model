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

nc = Dataset('Summit_basecase_Javis.nc','r')
# nc = Dataset('Summit_PD_fexp_0.35.nc','r')
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
'''
fig,(ax,ax0) = plt.subplots(nrows=2,sharex=True,figsize=(8,8))
plt.subplots_adjust(hspace=0)
plt.subplots_adjust(wspace=0)

#ax.set_ylim([-5,20])
ax.set_xlim([0,190])
ax.tick_params(direction="in")
eps_surface[eps_surface==0] = np.nan
ax.plot(np.arange(len(FD_d15N)),eps_surface,color='k',label=r'$^{15}$$\varepsilon$${_p}$')
ax.set_ylabel(r'10$^{3} \times$ $\varepsilon$${_p}$',fontsize=14,color='k')
ax.text(3,-60,'(a)',fontsize=14)
ax1 = ax.twinx()
ax1.plot(np.arange(len(FD_d15N)),(1-frac_rem)*100,'r--',label=r'f$_{loss}$')
ax1.set_ylabel('loss fraction (%)',fontsize=14,color='r')
lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax1.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2,loc=(0.47,0.75),fontsize=12,ncol=1,frameon=False)
ax1.tick_params(axis='y', colors='r')
ax.set_ylim([-90,-55])
for i in range(4):
    ax.fill([8+52*i,21+52*i,21+52*i,8+52*i],[-90,-90,-55,-55],'g',alpha=0.2)
# summer
for i in range(4):
    ax.fill([21+52*i,33+52*i,33+52*i,21+52*i],[-90,-90,-55,-55],'r',alpha=0.2)
# autumn
for i in range(3):
    ax.fill([33+52*i,47+52*i,47+52*i,33+52*i],[-90,-90,-55,-55],'b',alpha=0.2)
#second subplot
# plot average PIF
# index of season bin
idx = [0,9,22,35,48,61,74,87,100,113,126,139,152,165,178,187]
PIF_ave = np.zeros((len(idx)-1))
PIF_std = np.zeros((len(idx)-1))
for i in range(len(PIF_ave)):
    PIF_ave[i] = PIF[idx[i]:idx[i+1]].mean()
    PIF_std [i]= PIF[idx[i]:idx[i+1]].std()
x = (np.array(idx[1:])+(np.array(idx[:-1])))/2
# ax0.errorbar(x,PIF_ave,yerr=PIF_std,fmt='o',ecolor='k',elinewidth=1,ms=5,mfc='wheat',mec='salmon',capsize=3,label='PIE_model') 
# data from Javis 2008
xob = [78,  91.5, 104.5, 117.5, 130.5, 143.5, 156.5, 169.5]
yob = [6.7, -3.13, 6.1, 10.22, 6.88, 1.46, 5.49, 18.37]
ax0.scatter(xob,yob,s=30,c='k',marker='*',label='PIE_ob')
ax0.plot(np.arange(len(FD_d15N)),PIF,color='k',label='PIE')
ax0.xaxis.set_major_locator(ticker.FixedLocator((np.arange(8)*26)))
ax0.xaxis.set_major_formatter(ticker.FixedFormatter(([0,26,0,26,0,26,0,26])))
ax0.set_xlabel('week of year',fontsize=14)
ax0.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N',fontsize=14)
ax0.legend(loc='upper center',frameon=False)
ax0.set_ylim([-5,20])
for i in range(4):
    ax0.fill([8+52*i,21+52*i,21+52*i,8+52*i],[-5,-5,20,20],'g',alpha=0.2)
# summer
for i in range(4):
    ax0.fill([21+52*i,33+52*i,33+52*i,21+52*i],[-5,-5,20,20],'r',alpha=0.2)
# autumn
for i in range(3):
    ax0.fill([33+52*i,47+52*i,47+52*i,33+52*i],[-5,-5,20,20],'b',alpha=0.2)

ax2 = ax.twiny()
ax2.set_xlim([0,190])
ax2.xaxis.set_major_locator(ticker.FixedLocator((np.arange(8)*26)))
ax2.xaxis.set_major_formatter(ticker.FixedFormatter((['','2004','','2005','','2006','','2007'])))

labels = [-5,0,5,10,15,'']
ax0.set_yticklabels(labels)
ax0.text(3,17,'(b)',fontsize=14)

plt.savefig('Figure3_used.png',dpi=300)
# plt.savefig("Figure 3.pdf", bbox_inches='tight')
'''
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 04 04:15:04 2019

@author: admin
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

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

fexp_str = ['0','025','05','075','1']
linsty = ['b-','m-','c-','g-','y-']
fexp = [0,0.25,0.5,0.75,1]
path = 'Summit_fexp_'

#observation data
with open('observation.txt','r') as f:
    data = np.genfromtxt(f,skip_header=1)
####################################################################
file_name = 'summit_model.png'
fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(9, 5))
plt.subplots_adjust(hspace=0.1)
plt.subplots_adjust(wspace=0.1)

#concentration
ax1.set_xticks([])
ax1.set_ylabel('Deposited month', fontsize=14, color='k')
ax1.set_yticks([-1.983,-1.564,-1.29,-0.917,-0.643,-0.269,0])
ax1.set_yticklabels(['Jul','Jan','Jul','Jan','Jul','Jan','Jul'],fontsize=14)

ax11 = ax1.twiny()
ax11.text(40,-0.07,'(a)',fontsize=14)
#ax11.plot(SN_conc[:2182], -depth[:2182], color='silver',linestyle='-',linewidth=1,alpha=0.5)
# plot the 3 cm interval depth profile
#ax11.plot(res[:,1], -res[:,0], color='k',linestyle='-',linewidth=1,markersize=4,label='model simulation')
ax11.set_xlabel('[NO$_{3}$$^{-}$], ng.g$^{-1}$', fontsize=14, color='k', ha="center", va="center", labelpad=14)
ax11.plot(data[:,1],-data[:,0]/100,'r*--',linewidth=1,label='observation')
#ax11.set_xticks([0,100,200,300,400])

#d15N
ax2.set_xticks([])
ax22 = ax2.twiny()
#ax22.plot(SN_d15N[40:2182], -depth[40:2182], 'k-',linewidth=1,label='model simulation')
ax22.set_xlabel(r'10$^{3} \times$ $\delta^{15}$N', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel
ax22.plot(data[:,2],-data[:,0]/100,'r*--',linewidth=1,label='observation')
ax22.text(-13,-0.07,'(b)',fontsize=14)
#ax22.set_xticks([-12,-8,-4,0,4,8,12])
#D17O
ax3.set_xticks([])
ax33 = ax3.twiny()
#ax33.plot(SN_D17O[40:2182], -depth[40:2182], 'k-',linewidth=1,label='model simulation')
ax33.set_xlabel(r'10$^{3} \times$ $\Delta^{17}$O', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel
ax33.plot(data[:,3],-data[:,0]/100,'r*--',linewidth=1,label='observation')
ax33.text(26.2,-0.07,'(c)',fontsize=14)

ax333 = ax3.twinx()
ax333.set_ylim(ax33.get_ylim())
#ax333.set_yticks(np.arange(-2.5,0.5,0.5))
ax333.xaxis.label.set_color('red')
ax333.set_ylabel('Depth (m)',fontsize=14)

# polt concentration
for i in range(5):
    file_path = path+fexp_str[i]+'.nc'
    nc = Dataset(file_path,'r')
    depth_profile = nc['outputs1']
    depth = depth_profile['SN_depth'][:]
    SN_conc = depth_profile['SN_conc'][:]
    SN_d15N = depth_profile['SN_d15N'][:]
    SN_D17O = depth_profile['SN_D17O'][:]
    SN_data = depth_profile['SN_data'][:]
    SN_thickness = depth_profile['SN_thickness'][:]
    nc.close()
    interval = 0.03
    res = cal_mean(interval,depth,SN_conc,SN_d15N,SN_D17O,end=2.1) # varied Fpri, photolysis
    ax11.plot(res[:,1], -res[:,0],linsty[i],linewidth=1,alpha=0.5)
    ax22.plot(SN_d15N[:2182], -depth[:2182],linsty[i],linewidth=1,label='f$_{exp}$ = %.2f'%fexp[i],alpha=0.5)
    ax33.plot(SN_D17O[:2182], -depth[:2182],linsty[i],linewidth=1,alpha=0.5)

ax22.legend(loc='best', bbox_to_anchor=(-0.35, -0.2, 0.5, 0.5),ncol=3,fontsize=10)

nc = Dataset('Summit_fexp_0.nc','r')
#nc = Dataset('Summit_basecase.nc','r')
# nc = Dataset('Summit_temperature_dependent_phi.nc','r')
depth_profile = nc['outputs1']
depth = depth_profile['SN_depth'][:]
SN_conc = depth_profile['SN_conc'][:]
nc.close()
ax11.plot(SN_conc[:2182], -depth[:2182], color='silver',linestyle='-',linewidth=1,alpha=0.5)

plt.savefig('snow_profile.png', bbox_inches='tight', dpi=600)

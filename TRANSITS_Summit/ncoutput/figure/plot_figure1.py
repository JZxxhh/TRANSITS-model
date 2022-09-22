# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 21:29:10 2020

@author: admin
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
# from mpl_toolkits.basemap import Basemap

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


#observation data
with open('observation.txt','r') as f:
    data = np.genfromtxt(f,skip_header=1)
    
with open('Fpri_partition.txt','r') as f:
    Fpri = np.genfromtxt(f,skip_header=0)

Fpri = Fpri*6.6e-6*52

####################################################################
file_name = 'snow_profile.eps'
fig, (ax1,ax2,ax3) = plt.subplots(1,3,sharey=True,figsize=(9, 5))
plt.subplots_adjust(hspace=0.1)
plt.subplots_adjust(wspace=0.1)

#concentration
# ax1.set_xticks([])
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
ax22.plot(data[:,2],-data[:,0]/100,'r*--',linewidth=1)
ax22.text(-13,-0.07,'(b)',fontsize=14)

# plot the observed surface snow nitrate d15N
# depth_sur = [-0.1345,-0.269,-0.456,-0.643,-0.78,-0.917,-1.1035,-1.29]
# d15N_sur = [-11.7,-13,-10.2,-2.8,-6.7,-10.9,-0.1,1.3]
# ax22.scatter(d15N_sur,depth_sur,s=20,c='k',marker='o')


#ax22.set_xticks([-12,-8,-4,0,4,8,12])
#D17O
ax3.set_xticks([])
ax33 = ax3.twiny()
#ax33.plot(SN_D17O[40:2182], -depth[40:2182], 'k-',linewidth=1,label='model simulation')
ax33.set_xlabel(r'10$^{3} \times$ $\Delta^{17}$O', fontsize=14, color='k', ha="center", va="center", labelpad=14) #set xlabel and ylabel
ax33.plot(data[:,3],-data[:,0]/100,'r*--',linewidth=1)
ax33.text(26.2,-0.07,'(c)',fontsize=14)

ax333 = ax3.twinx()
ax333.set_ylim(ax33.get_ylim())
#ax333.set_yticks(np.arange(-2.5,0.5,0.5))
ax333.xaxis.label.set_color('red')
ax333.set_ylabel('Depth (m)',fontsize=14)

##########
#plot modelled profile
# file_path = 'Summit_fexp_015.nc'
# file_path = 'Summit_varied_Fpri_no_photolysis.nc'
file_path = 'Summit_PD_fexp_0.35.nc'
# file_path = 'Summit_varied_Fpri.nc'
nc = Dataset(file_path,'r')
depth_profile = nc['outputs1']
depth = depth_profile['SN_depth'][:]
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_D17O = depth_profile['SN_D17O'][:]
SN_data = depth_profile['SN_data'][:]
SN_thickness = depth_profile['SN_thickness'][:]
nc.close()

# calculate the Fpri flux
Fpri_flux = np.zeros_like(SN_conc)
for i in range(len(Fpri_flux)):
    Fpri_flux[i] = Fpri[SN_data[i]%52]

file_path = 'Summit_varied_Fpri_no_photolysis_01.nc'
nc = Dataset(file_path,'r')
depth_profile = nc['outputs1']
SN_d15N1 = depth_profile['SN_d15N'][:]

interval = 0.03
res = cal_mean(interval,depth,SN_conc,SN_d15N,SN_D17O,end=2.1) # varied Fpri, photolysis
ax11.plot(res[:,1], -res[:,0],color='b',linestyle='-',
          linewidth=1,label='model (varied $\delta^{15} $N_F$_{pri}$)')
ax22.plot(SN_d15N[20:2182], -depth[20:2182],color='b',linestyle='-',linewidth=1)
# ax22.plot(SN_d15N1[20:2182], -depth[20:2182],'g:',linewidth=2,label=r'$\delta^{15}$N_F$_{pri}$')
ax33.plot(SN_D17O[:2182], -depth[:2182],'b-',linewidth=1)
ax11.plot(SN_conc[:2182], -depth[:2182], color='silver',linestyle='-',linewidth=1,alpha=0.5)

# ax11.legend(loc='lower right',fontsize=10)

# ax22.plot(SN_d15N[20:2182]-SN_d15N1[20:2182], -depth[20:2182],
#           color='k',linestyle=(0, (5, 5)),linewidth=1,label='model (constant $\delta^{15}$N_F$_{pri}$)')

ax1.plot(Fpri_flux[:2182]*1e6,-depth[:2182],'k-',linewidth=1,label=r'F$_{pri}$')
ax1.set_xlabel(r'$kgN.m^{-2}.yr^{-1}$', fontsize=14, color='k')
ax1.text(8.1,-2.3,r'$\times10^{-6}$',fontsize=8)

lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
order = [0,3,2,4,1]
lines, labels = [lines[idx] for idx in order],[labels[idx] for idx in order]
fig.legend(lines, labels, loc = (0.36, 0.01), ncol=3)

# plt.savefig(file_name,format='eps')

# plt.savefig("Figure 1.pdf", bbox_inches='tight')













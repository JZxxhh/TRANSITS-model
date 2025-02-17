# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 10:25:50 2024

@author: pb130
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt



dt = 7*24*3600 
# nc = Dataset('test_SE_Dome.nc','r')
nc = Dataset('present_day_without_ozone_hole.nc','r')

fout = nc['parameters']['fout'][:]
weekly_profile = nc['outputs2']

# weekly variables
FP_mass = weekly_profile['FP_mass'][:]*dt
FP_d15N = weekly_profile['FP_d15N'][:]
FP_D17O = weekly_profile['FP_D17O'][:]
FD_mass = weekly_profile['FD_mass'][:]*dt
FD_d15N = weekly_profile['FD_d15N'][:]
FD_D17O = weekly_profile['FD_D17O'][:]
eps_surface = weekly_profile['eps_surface'][:]
eps_surface = weekly_profile['eps_surface'][:]
frac_rem = weekly_profile['frac_rem'][:]
surface_conc = weekly_profile['surface_conc'][:-25]
surface_d15N = weekly_profile['surface_d15N'][:-25]
surface_D17O = weekly_profile['surface_D17O'][:-25]

depth_profile = nc['outputs1']
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_data = depth_profile['SN_data'][:]
SN_D17O = depth_profile['SN_D17O'][:]
SN_thickness = depth_profile['SN_thickness'][:]
nc.close()

depth = SN_thickness.cumsum()

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
    FA_D17O[i] = (FA[idx]*SN_D17O[idx]).sum()/FA[idx].sum() # note D17O is inaccurate here becasue I didn't modify the oxidants in model.

# using the annual value in 5th year

# idx = 260: 312
FA_annual = FA_mass[260: 312].mean()*52*dt
FA_d15N_annual = (FA_mass[260: 312]*FA_d15N[260: 312]).sum()/FA_mass[260: 312].sum()

print('%.2f  %.2f'%(FA_d15N_annual,(1-FA_annual/Fpri)*100))

SN_ob = np.genfromtxt('WD_snowpack_nitrate.txt',skip_header=1)
# comapring with the snowpack observation

fig, (ax0,ax1) = plt.subplots(1,2,figsize=(6, 5))
plt.subplots_adjust(wspace=0.1)

ax01 = ax0.twiny()
ax01.plot(SN_conc,-depth,'k-',linewidth=2,label='Model')
ax01.set_ylim([-2,0])
ax01.set_xlim([0,300])

# Now load the snow nitrate concentration profile observed from 
# a thick, 2m depth snowpack and 3 shallow firns from Masclin 2013 ACP paper
SN_ob1 = np.genfromtxt('WD_snowpack_nitrate.txt',skip_header=1) 

SN_ob2 = np.genfromtxt('WD_snow_nitrate_shallow_firn.txt',skip_header=1) # note the unit of concentration is in miuM


ax01.scatter(SN_ob1[:,1],-SN_ob1[:,0]/100,s=28,c='b',marker='*',label='Kreutz et al. (2008)')

ax01.scatter(SN_ob2[:,0]*62,SN_ob2[:,1]/100,s=20,c='r',
             marker='^',label='Masclin et al. (2013)')
ax01.scatter(SN_ob2[:,2]*62,SN_ob2[:,3]/100,s=20,c='salmon',
             marker='v')
ax01.scatter(SN_ob2[:,4]*62,SN_ob2[:,5]/100,s=20,c='plum',
             marker='>')

ax01.legend(loc=(0.2,0.25),fontsize=9)

ax01.set_xlabel('[NO$_{3}$$^{-}$], ng.g$^{-1}$',
                fontsize=12, color='k', ha="center",
                va="center", labelpad=14)

ax0.set_ylabel('Depth (m)',fontsize=12)
ax0.set_xticks([])

ax01.fill([0,0,300,300],[0,-0.8,-0.8,0],'yellow',alpha=0.3)

ax01.set_xticks([0,100,200])

# Panel 2: loss fraction and 
# weekly depth
dh = 200/52/SN_d # accumulated thickness per week
depth_weeky = (np.ones((nb_tt,)).cumsum()-0.5)*dh # mid depth


ax11 = ax1.twiny()

ax11.plot((1-frac_rem[::-1])*100,-depth_weeky ,'r-',linewidth=2,label='Lost fraction')


ax12 = ax1.twinx()
ax12.plot(FA_d15N[::-1],-depth_weeky ,'b--',linewidth=2,label='Model')


ax11.set_ylim([-2,0])
ax12.set_ylim([-2,0])
ax11.set_yticks([])


ax11.set_xlabel('Nitrate lost fraction (%)',
                fontsize=12, color='r', ha="center",
                va="center", labelpad=14)

ax1.set_xlabel(r'10$^{3} \times$ $\delta^{15}$N(NO$_{3}$$^{-}$)', 
                fontsize=12, color='b', ha="center",
                va="center", labelpad=14)

ax11.set_xlim([0,80])
ax11.fill([0,0,80,80],[0,-0.8,-0.8,0],'yellow',alpha=0.3)


ax0.text(-0.1,1.1,'(a)',fontsize=14,transform=ax0.transAxes)
ax1.text(-0.1,1.1,'(b)',fontsize=14,transform=ax1.transAxes)

plt.savefig('WD_snowpack_profile_without_ozonehole.png',dpi=300,bbox_inches='tight')





























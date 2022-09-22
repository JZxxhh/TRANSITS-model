# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:40:23 2020

@author: admin
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


dt = 7*24*3600 
nc = Dataset('Summit_varied_Fpri_d15N.nc','r')

weekly_profile = nc['outputs2']
# weekly variables
FP_mass = weekly_profile['FP_mass'][:-25]*dt
FP_d15N = weekly_profile['FP_d15N'][:-25]
FP_D17O = weekly_profile['FP_D17O'][:-25]
FD_mass = weekly_profile['FD_mass'][:-25]*dt
FD_d15N = weekly_profile['FD_d15N'][:-25]
FD_D17O = weekly_profile['FD_D17O'][:-25]
eps_surface = weekly_profile['eps_surface'][:-25]
eps_surface = weekly_profile['eps_surface'][:-25]
frac_rem = weekly_profile['frac_rem'][:-25]
#surface_conc = weekly_profile['surface_conc'][:-25]
#surface_d15N = weekly_profile['surface_d15N'][:-25]
#surface_D17O = weekly_profile['surface_D17O'][:-25]

depth_profile = nc['outputs1']
SN_conc = depth_profile['SN_conc'][:]
SN_d15N = depth_profile['SN_d15N'][:]
SN_data = depth_profile['SN_data'][:]
SN_D17O = depth_profile['SN_D17O'][:]
SN_thickness = depth_profile['SN_thickness'][:]
nc.close()

# find the weekly preserved conc, d15N and D17O
res = np.zeros((52,4))
for week in range(52):
    conc = []
    d15N = []
    D17O = []
    thickness = []
    for layer in range(len(SN_data)):
        if SN_data[layer] == week+52:
            conc.append(SN_conc[layer])
            d15N.append(SN_d15N[layer])
            D17O.append(SN_D17O[layer])
            thickness.append(SN_thickness[layer])
    res[week,0] = (np.array(conc)*np.array(thickness)).sum()/np.array(thickness).sum()
    res[week,1] = (np.array(conc)*np.array(thickness)*np.array(d15N)).sum()/(np.array(conc)*np.array(thickness)).sum()
    res[week,2] = (np.array(conc)*np.array(thickness)*np.array(D17O)).sum()/(np.array(conc)*np.array(thickness)).sum()
    res[week,3] = np.array(thickness).sum() #accumulated thickness

# calculate the average conc, d15N and D17O of a intergral year
# week 52:104 is used
for i in range(len(SN_data)):
    if SN_data[i] == 103 and SN_data[i-1] !=103:
        idx1 = i
    if SN_data[i+1] !=52 and SN_data[i] ==52:
        idx2 = i
        break
conca = (SN_conc[idx1:idx2+1]*SN_thickness[idx1:idx2+1]).sum()/SN_thickness[idx1:idx2+1].sum()
d15Na = (SN_d15N[idx1:idx2+1]*SN_conc[idx1:idx2+1]*SN_thickness[idx1:idx2+1]).sum()/(SN_conc[idx1:idx2+1]*SN_thickness[idx1:idx2+1]).sum()
D17Oa = (SN_D17O[idx1:idx2+1]*SN_conc[idx1:idx2+1]*SN_thickness[idx1:idx2+1]).sum()/(SN_conc[idx1:idx2+1]*SN_thickness[idx1:idx2+1]).sum()
# calculate FA
FA = (SN_conc[idx1:idx2+1]*SN_thickness[idx1:idx2+1]*380*1000*1e-12*14/62).sum()



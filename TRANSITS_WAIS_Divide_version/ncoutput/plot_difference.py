# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 12:28:52 2024

@author: pb130
"""

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt


# observation results
yo = [6.3, 32.9, 26.6] # observations of Holocene, glacial average and their difference
err_ob = [2.5, 5.1, 5.7] # standard error of observations

# model results from statistic.py
# Holocene archived d15N and loss fraction: 9.33  12.72
# glacial archived d15N and loss fraction: 31.08  31.81
# glacial archived d15N and loss fraction with phi(T): 26.30  28.16
xob = [0.7, 2.3, 4.5]
yob = yo
err_yob = err_ob

xmd = [1.3, 2.9, 5.1] # location of modeled value
ymd = [9.33, 31.08, 21.75] # modeled Holocene, glacial d15N and their difference
floss_md = ['(12.7)','(31.8)','']

xmd_T = [3.5, 5.7]
ymd_T = [26.30, 16.97] # modeled glacial d15N and their difference with considering temperature effect on phi
floss_md_T = ['(28.2)','']

bar_width = 0.6
error_param = dict(elinewidth=3, ecolor='k', capsize=4)

fig, ax = plt.subplots(figsize=(8,6))

bar1 = ax.bar(xob, yob, bar_width, color='r',label='Observation',
              yerr=err_yob,error_kw=error_param,alpha=0.5)
bar2 = ax.bar(xmd, ymd, bar_width, color='b',label=r'Model (constant $\phi$)',alpha=0.5)
bar3 = ax.bar(xmd_T, ymd_T, bar_width, color='y',label=r'Model (with $\phi$(T))',alpha=0.5)

# ax.bar_label(bar2, labels=floss_md, padding=2,fontsize=11)
# ax.bar_label(bar3, labels=floss_md_T, padding=2,fontsize=11)

ax.legend(loc='upper right',fontsize=12)

ax.set_xticks([1, 2.9, 5.1],['Holocene','Glacial','Glacial-Holocene\ndifference'],fontsize=14)
ax.set_ylim([0,45])

ax.set_ylabel(r'10$^{3} \times$ $\delta^{15}$N(NO$_{3}$$^{-}$)',fontsize=14)
plt.savefig('Glacial_interglacial_difference_d15N.png',dpi=300,bbox_inches='tight')


























